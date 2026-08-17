#!/bin/bash

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
# FORCE MAC NATIVE JAVA 21 PATH (Homebrew)
# ==========================================
JAVA_21_EXEC="/opt/homebrew/opt/openjdk@21/bin/java"

if [ ! -f "$JAVA_21_EXEC" ]; then
    echo "❌ Error: Java 21 was not found at $JAVA_21_EXEC"
    exit 1
fi

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

BATCH_FILE_NAME="mac_bulk_temp.bat"
ABS_BATCH_FILE="$SCRIPT_DIR/$BATCH_FILE_NAME"

# ==========================================
# PROCESSING LOOP
# ==========================================
for folder in "$PARENT_DIR"/*/; do
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"
        
        # 1. FIX: Find only FASTA files matching pPSBA_[0-9].fa exactly (ignores suffixes)
        # Using Mac-compatible find syntax with extended regular expressions
        fa_file=$(find -E "$folder" -maxdepth 1 -type f -regex '.*/pPSBA_[0-9]+\.fa' | head -n 1)
        
        if [ -z "$fa_file" ]; then
            echo "   ⚠️ Skipping: No valid pPSBA_[NUMBER].fa reference file found in $folder"
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

        # 4. Create an output directory for snapshots inside this folder
        out_dir="${folder}Snapshots"
        mkdir -p "$out_dir"

        # 5. Build individual batch commands for this folder's contents
        echo "new" > "$ABS_BATCH_FILE"
        echo "genome $fa_file" >> "$ABS_BATCH_FILE"
        echo "snapshotDirectory $out_dir" >> "$ABS_BATCH_FILE"
        echo "maxPanelHeight 1000" >> "$ABS_BATCH_FILE"

        bam_count=0
        for bam_file in "$bam_folder"/*.bam; do
            if [ -f "$bam_file" ]; then
                filename=$(basename "$bam_file" .bam)
                
                echo "load $bam_file" >> "$ABS_BATCH_FILE"
                echo "goto $SEQUENCE_ID" >> "$ABS_BATCH_FILE"
                echo "snapshot ${filename}_snapshot.png" >> "$ABS_BATCH_FILE"
                echo "remove $bam_file" >> "$ABS_BATCH_FILE"
                ((bam_count++))
            fi
        done

        echo "exit" >> "$ABS_BATCH_FILE"

        # 6. Run IGV using Classpath strings
        if [ $bam_count -gt 0 ]; then
            echo "   🎯 Valid Reference Identified: $(basename "$fa_file")"
            echo "   📸 Sequence Targeted: [$SEQUENCE_ID]"
            echo "   📸 Generating $bam_count snapshots using IGV..."
            
            cd "$IGV_DIR" || exit
            "$JAVA_21_EXEC" -Xmx4000m \
                -Djava.net.useSystemProxies=true \
                -Dapple.laf.useScreenMenuBar=true \
                -cp "lib/*" \
                org.broad.igv.ui.Main --batch "$ABS_BATCH_FILE"
            cd "$SCRIPT_DIR" || exit
                
        else
            echo "   ⚠️ No .bam files found inside $bam_folder"
        fi
        
        # Clean up temporary batch file
        rm -f "$ABS_BATCH_FILE"
    fi
done

echo "🎉 All directories completed!"
