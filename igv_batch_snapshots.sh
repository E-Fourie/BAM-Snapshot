#!/bin/bash
#
# igv_batch_snapshots.sh
# Mac / Linux batch IGV snapshot generator.
# For native Windows, use igv_batch_snapshots.ps1 instead.

# ==========================================
# ROBUST IGV DIRECTORY LOOKUP
# ==========================================
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
IGV_DIR=""

for d in "$SCRIPT_DIR"/[Ii][Gg][Vv]*/; do
    if [ -d "${d}lib" ]; then
        IGV_DIR="${d%/}"
        break
    fi
done

if [ -z "$IGV_DIR" ] && [ -d "$SCRIPT_DIR/lib" ]; then
    IGV_DIR="$SCRIPT_DIR"
fi

if [ -z "$IGV_DIR" ]; then
    echo "❌ Error: Could not find the IGV folder containing a 'lib' directory inside $SCRIPT_DIR"
    exit 1
fi

echo "🔍 Found IGV directory at: $IGV_DIR"

# ==========================================
# OS DETECTION (drives java lookup + find syntax)
# ==========================================
OS_NAME=$(uname -s)

case "$OS_NAME" in
    Darwin) PLATFORM="mac" ;;
    Linux)  PLATFORM="linux" ;;
    *)
        echo "❌ Error: Unsupported platform '$OS_NAME'. Use igv_batch_snapshots.ps1 on native Windows."
        exit 1
        ;;
esac

# ==========================================
# JAVA LOOKUP
# Preference order:
#   1. A JDK bundled alongside IGV itself (jdk*/bin/java under IGV_DIR).
#   2. Common Homebrew Java 21 install locations (mac).
#   3. 'java' already on the PATH.
# ==========================================
JAVA_EXEC=""

bundled_java=$(find "$IGV_DIR" -maxdepth 3 -type f -iname "java" -path "*/bin/*" 2>/dev/null | head -n 1)
if [ -n "$bundled_java" ] && [ -x "$bundled_java" ]; then
    JAVA_EXEC="$bundled_java"
fi

if [ -z "$JAVA_EXEC" ] && [ "$PLATFORM" = "mac" ]; then
    for candidate in \
        "/opt/homebrew/opt/openjdk@21/bin/java" \
        "/usr/local/opt/openjdk@21/bin/java" \
        "/opt/homebrew/opt/openjdk/bin/java" \
        "/usr/local/opt/openjdk/bin/java"
    do
        if [ -x "$candidate" ]; then
            JAVA_EXEC="$candidate"
            break
        fi
    done
fi

if [ -z "$JAVA_EXEC" ] && command -v java >/dev/null 2>&1; then
    JAVA_EXEC=$(command -v java)
fi

if [ -z "$JAVA_EXEC" ]; then
    echo "❌ Error: Could not find a Java executable (checked bundled JDK, Homebrew, and PATH)."
    echo "   Install Java 21 (e.g. 'brew install openjdk@21' on Mac, or your package manager on Linux)."
    exit 1
fi

echo "☕ Using Java at: $JAVA_EXEC"

# ==========================================
# COMMAND-LINE ARGUMENT CHECK
# ==========================================
if [ -z "$1" ]; then
    echo "❌ Error: No directory path provided."
    echo "Usage: $0 /path/to/your/main/directory"
    exit 1
fi

PARENT_DIR=$(cd "$1" 2>/dev/null && pwd)

if [ -z "$PARENT_DIR" ] || [ ! -d "$PARENT_DIR" ]; then
    echo "❌ Error: The directory '$1' does not exist."
    exit 1
fi

echo "🚀 Starting batch IGV snapshots..."
echo "📂 Processing Main Directory: $PARENT_DIR"
echo "------------------------------------------------"

BATCH_FILE_NAME="bulk_temp.bat"
ABS_BATCH_FILE="$SCRIPT_DIR/$BATCH_FILE_NAME"

# ==========================================
# FASTA LOOKUP (platform-specific regex syntax)
#   Matches any .fa / .fasta file, excluding anything containing "_clean".
# ==========================================
find_reference_fasta() {
    local folder="$1"
    if [ "$PLATFORM" = "mac" ]; then
        find -E "$folder" -maxdepth 1 -type f -iregex '.*\.(fa|fasta)' ! -iname '*_clean*' | head -n 1
    else
        find "$folder" -maxdepth 1 -type f -regextype posix-extended -iregex '.*\.(fa|fasta)' ! -iname '*_clean*' | head -n 1
    fi
}

# ==========================================
# PROCESSING LOOP
# ==========================================
for folder in "$PARENT_DIR"/*/; do
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"

        # 1. Match ANY .fa/.fasta reference file, skipping *_clean* variants.
        fa_file=$(find_reference_fasta "$folder")

        if [ -z "$fa_file" ]; then
            echo "   ⚠️ Skipping: No valid (non-_clean) .fa/.fasta reference file found in $folder"
            continue
        fi

        # 2. Check if the subfolder "BAM" exists
        bam_folder="${folder}BAM"
        if [ ! -d "$bam_folder" ]; then
            echo "   ⚠️ Skipping: No 'BAM' subfolder found in $folder"
            continue
        fi

        # 3. Dynamically extract the sequence ID header name from the valid .fa file
        SEQUENCE_ID=$(grep "^>" "$fa_file" | head -n 1 | awk '{print $1}' | sed 's/>//')

        if [ -z "$SEQUENCE_ID" ]; then
            echo "   ⚠️ Skipping: Could not read a valid sequence ID header from $fa_file"
            continue
        fi

        # 4. Real output directory for the final snapshots
        out_dir="${folder}Snapshots"
        mkdir -p "$out_dir"

        # ==========================================
        # FIX: IGV's batch-file parser has a long-standing, unfixed bug —
        # it does NOT strip quote characters from quoted paths (see
        # igvteam/igv#1258), so quoting a path with spaces just makes IGV
        # look for a filename that literally includes the quote marks and
        # fail. That failed "genome" load is what leaves
        # GenomeManager.getCurrentGenome() null and crashes the very next
        # command. There is no working quoting syntax for this — spaces in
        # a path are simply unsupported by IGV batch mode.
        #
        # The reliable workaround is to copy the fasta/BAM/index files IGV
        # needs into a throwaway, guaranteed-space-free temp directory, and
        # reference *those* paths (completely unquoted) in the batch file.
        # ==========================================
        WORK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/igv_batch.XXXXXX")

        fa_basename=$(basename "$fa_file")
        cp "$fa_file" "$WORK_DIR/$fa_basename"
        chmod u+w "$WORK_DIR/$fa_basename"
        if [ -f "${fa_file}.fai" ]; then
            cp "${fa_file}.fai" "$WORK_DIR/${fa_basename}.fai"
            chmod u+w "$WORK_DIR/${fa_basename}.fai"
        elif command -v samtools >/dev/null 2>&1; then
            samtools faidx "$WORK_DIR/$fa_basename" 2>/dev/null
        fi

        work_bam_dir="$WORK_DIR/BAM"
        mkdir -p "$work_bam_dir"
        work_snap_dir="$WORK_DIR/Snapshots"
        mkdir -p "$work_snap_dir"

        # 5. Build batch commands using the copied, space-free paths.
        #    No quotes needed (and none used — see note above).
        # NOTE: "genome" must come before "new". IGV's "new" command calls
        # newSession(), which resets the view using whatever genome is
        # currently active. On a machine with no cached "last used genome"
        # (e.g. a fresh IGV preferences directory), nothing is loaded yet
        # at batch-script start, so "new" throws a NullPointerException
        # (GenomeManager.getCurrentGenome() is null) before the script's
        # own "genome" line ever runs. Loading the genome first avoids this
        # regardless of what IGV had cached from a previous run.
        echo "genome $WORK_DIR/$fa_basename" > "$ABS_BATCH_FILE"
        # IGV loads the genome asynchronously in batch mode; without a
        # pause the next command can run before the genome finishes
        # loading. This gives it time.
        echo "setSleepInterval 2000" >> "$ABS_BATCH_FILE"
        echo "new" >> "$ABS_BATCH_FILE"
        echo "snapshotDirectory $work_snap_dir" >> "$ABS_BATCH_FILE"
        echo "maxPanelHeight 1000" >> "$ABS_BATCH_FILE"

        bam_count=0
        for bam_file in "$bam_folder"/*.bam; do
            if [ -f "$bam_file" ]; then
                bam_basename=$(basename "$bam_file")
                cp "$bam_file" "$work_bam_dir/$bam_basename"
                chmod u+w "$work_bam_dir/$bam_basename"
                if [ -f "${bam_file}.bai" ]; then
                    cp "${bam_file}.bai" "$work_bam_dir/${bam_basename}.bai"
                    chmod u+w "$work_bam_dir/${bam_basename}.bai"
                fi

                echo "load $work_bam_dir/$bam_basename" >> "$ABS_BATCH_FILE"
                echo "goto $SEQUENCE_ID" >> "$ABS_BATCH_FILE"
                echo "snapshot ${bam_basename%.bam}_snapshot.png" >> "$ABS_BATCH_FILE"
                echo "remove $work_bam_dir/$bam_basename" >> "$ABS_BATCH_FILE"
                ((bam_count++))
            fi
        done

        echo "exit" >> "$ABS_BATCH_FILE"

        # 6. Run IGV
        if [ $bam_count -gt 0 ]; then
            echo "   🎯 Valid Reference Identified: $fa_basename"
            echo "   📸 Sequence Targeted: [$SEQUENCE_ID]"
            echo "   📸 Generating $bam_count snapshots using IGV..."

            cd "$IGV_DIR" || exit
            MAC_FLAGS=()
            if [ "$PLATFORM" = "mac" ]; then
                MAC_FLAGS=(-Dapple.laf.useScreenMenuBar=true)
            fi
            "$JAVA_EXEC" -Xmx4000m \
                -Djava.net.useSystemProxies=true \
                "${MAC_FLAGS[@]}" \
                -cp "lib/*" \
                org.broad.igv.ui.Main --batch "$ABS_BATCH_FILE"
            cd "$SCRIPT_DIR" || exit

            # 7. Copy the resulting snapshots back to the real per-sample folder.
            shopt -s nullglob
            snaps=("$work_snap_dir"/*.png)
            shopt -u nullglob
            if [ ${#snaps[@]} -gt 0 ]; then
                cp "${snaps[@]}" "$out_dir/"
                echo "   ✅ Copied ${#snaps[@]} snapshot(s) to $out_dir"
            else
                echo "   ⚠️ No snapshot images were produced — check igv.log in $IGV_DIR"
            fi
        else
            echo "   ⚠️ No .bam files found inside $bam_folder"
        fi

        # Clean up
        rm -f "$ABS_BATCH_FILE"
        rm -rf "$WORK_DIR"
    fi
done

echo "🎉 All directories completed!"
