<#
.SYNOPSIS
    igv_batch_snapshots.ps1 - Windows batch IGV snapshot generator.
    For Mac/Linux, use igv_batch_snapshots.sh instead.

.USAGE
    Right-click -> "Run with PowerShell", or from a PowerShell prompt:
        .\igv_batch_snapshots.ps1 "C:\path\to\your\main\directory"

    If Windows blocks the script (unsigned script policy), run once:
        Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
    then re-run the command above in the same window.
#>

param(
    [Parameter(Mandatory = $true, Position = 0)]
    [string]$ParentDirArg
)

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

function Fail($msg) {
    Write-Host "$msg" -ForegroundColor Red
    exit 1
}

# ==========================================
# ROBUST IGV DIRECTORY LOOKUP
# ==========================================
$IgvDir = $null
$candidateDirs = Get-ChildItem -Path $ScriptDir -Directory -Filter "IGV*" -ErrorAction SilentlyContinue
foreach ($d in $candidateDirs) {
    if (Test-Path (Join-Path $d.FullName "lib")) {
        $IgvDir = $d.FullName
        break
    }
}
if (-not $IgvDir -and (Test-Path (Join-Path $ScriptDir "lib"))) {
    $IgvDir = $ScriptDir
}
if (-not $IgvDir) {
    Fail "Error: Could not find the IGV folder containing a 'lib' directory inside $ScriptDir"
}
Write-Host "Found IGV directory at: $IgvDir"

# ==========================================
# JAVA LOOKUP
# Preference order:
#   1. A JDK bundled alongside IGV itself (jdk*\bin\java.exe under $IgvDir).
#   2. Common Windows Java 21 install locations (Adoptium/Temurin, Oracle).
#   3. 'java' already on the PATH.
# ==========================================
$JavaExec = $null

$bundledJava = Get-ChildItem -Path $IgvDir -Recurse -Depth 3 -Filter "java.exe" -ErrorAction SilentlyContinue |
    Where-Object { $_.DirectoryName -match '\\bin$' } |
    Select-Object -First 1
if ($bundledJava) {
    $JavaExec = $bundledJava.FullName
}

if (-not $JavaExec) {
    $winCandidates = @(
        "$Env:ProgramFiles\Eclipse Adoptium\jdk-21*\bin\java.exe",
        "$Env:ProgramFiles\Java\jdk-21*\bin\java.exe",
        "$Env:ProgramFiles\Microsoft\jdk-21*\bin\java.exe"
    )
    foreach ($pattern in $winCandidates) {
        $found = Get-ChildItem -Path $pattern -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($found) { $JavaExec = $found.FullName; break }
    }
}

if (-not $JavaExec) {
    $onPath = Get-Command java.exe -ErrorAction SilentlyContinue
    if ($onPath) { $JavaExec = $onPath.Source }
}

if (-not $JavaExec) {
    Fail "Error: Could not find java.exe (checked bundled JDK, common install paths, and PATH). Install Java 21 (e.g. Eclipse Temurin) and try again."
}
Write-Host "Using Java at: $JavaExec"

# ==========================================
# COMMAND-LINE ARGUMENT CHECK
# ==========================================
if (-not (Test-Path $ParentDirArg -PathType Container)) {
    Fail "Error: The directory '$ParentDirArg' does not exist."
}
$ParentDir = (Resolve-Path $ParentDirArg).Path

Write-Host "Starting batch IGV snapshots..."
Write-Host "Processing Main Directory: $ParentDir"
Write-Host "------------------------------------------------"

$BatchFile = Join-Path $ScriptDir "bulk_temp.bat"

# ==========================================
# A SPACE-FREE SCRATCH ROOT
# IGV's batch parser has a long-standing, unfixed bug: it does not strip
# quote characters from quoted paths (igvteam/igv#1258), so quoting a path
# with spaces just makes IGV look for a filename that literally includes
# the quote marks and fail to load — which is what leaves
# GenomeManager.getCurrentGenome() null and crashes the very next command.
# There's no working quoting syntax for this; spaces are simply unsupported.
# So every file IGV touches gets copied into a guaranteed space-free
# scratch folder first, and the batch file references those copies,
# completely unquoted. $Env:TEMP can itself contain spaces (e.g. a
# Windows username with a space in it), so prefer $Env:ProgramData instead.
# ==========================================
$ScratchRoot = Join-Path $Env:ProgramData "PSBA_IGV_scratch"
try {
    New-Item -ItemType Directory -Path $ScratchRoot -Force -ErrorAction Stop | Out-Null
} catch {
    $ScratchRoot = $Env:TEMP
}
if ($ScratchRoot -match '\s') {
    Write-Host "Warning: scratch path '$ScratchRoot' contains a space; IGV batch mode may still fail. Consider setting `$ScratchRoot in the script to a space-free folder such as C:\IGVScratch." -ForegroundColor Yellow
}

# ==========================================
# PROCESSING LOOP
# ==========================================
$folders = Get-ChildItem -Path $ParentDir -Directory
foreach ($folderItem in $folders) {
    $folder = $folderItem.FullName
    Write-Host "Processing folder: $folder"

    # 1. Match ANY .fa/.fasta reference file, skipping *_clean* variants.
    $faFile = Get-ChildItem -Path $folder -File -ErrorAction SilentlyContinue |
        Where-Object { $_.Extension -match '^\.(fa|fasta)$' -and $_.Name -notmatch '(?i)_clean' } |
        Select-Object -First 1

    if (-not $faFile) {
        Write-Host "   Skipping: No valid (non-_clean) .fa/.fasta reference file found in $folder" -ForegroundColor Yellow
        continue
    }

    # 2. Check if the subfolder "BAM" exists
    $bamFolder = Join-Path $folder "BAM"
    if (-not (Test-Path $bamFolder -PathType Container)) {
        Write-Host "   Skipping: No 'BAM' subfolder found in $folder" -ForegroundColor Yellow
        continue
    }

    # 3. Dynamically extract the sequence ID header name from the valid .fa file
    $headerLine = Get-Content $faFile.FullName -TotalCount 50 | Where-Object { $_ -match '^>' } | Select-Object -First 1
    if (-not $headerLine) {
        Write-Host "   Skipping: Could not read a valid sequence ID header from $($faFile.FullName)" -ForegroundColor Yellow
        continue
    }
    $SequenceId = ($headerLine.TrimStart('>') -split '\s+')[0]

    # 4. Real output directory for the final snapshots
    $outDir = Join-Path $folder "Snapshots"
    New-Item -ItemType Directory -Path $outDir -Force | Out-Null

    # 5. Set up a per-sample space-free scratch folder and copy in only
    #    what IGV needs: the reference fasta (+ .fai if present) and the
    #    BAM files (+ .bai if present).
    $workDir = Join-Path $ScratchRoot ([guid]::NewGuid().ToString())
    $workBamDir = Join-Path $workDir "BAM"
    $workSnapDir = Join-Path $workDir "Snapshots"
    New-Item -ItemType Directory -Path $workBamDir -Force | Out-Null
    New-Item -ItemType Directory -Path $workSnapDir -Force | Out-Null

    $workFasta = Join-Path $workDir $faFile.Name
    Copy-Item -Path $faFile.FullName -Destination $workFasta -Force
    Set-ItemProperty -Path $workFasta -Name IsReadOnly -Value $false -ErrorAction SilentlyContinue
    $faiSource = "$($faFile.FullName).fai"
    if (Test-Path $faiSource) {
        $workFai = "$workFasta.fai"
        Copy-Item -Path $faiSource -Destination $workFai -Force
        Set-ItemProperty -Path $workFai -Name IsReadOnly -Value $false -ErrorAction SilentlyContinue
    }

    $bamFiles = Get-ChildItem -Path $bamFolder -Filter "*.bam" -File -ErrorAction SilentlyContinue

    # 6. Build batch commands using the copied, space-free paths.
    #    No quotes needed (and none used — see note above on the IGV bug).
    # NOTE: "genome" must come before "new". IGV's "new" command calls
    # newSession(), which resets the view using whatever genome is
    # currently active. On a machine with no cached "last used genome"
    # (e.g. a fresh IGV preferences directory), nothing is loaded yet at
    # batch-script start, so "new" throws a NullPointerException
    # (GenomeManager.getCurrentGenome() is null) before the script's own
    # "genome" line ever runs. Loading the genome first avoids this
    # regardless of what IGV had cached from a previous run.
    $lines = @()
    $lines += "genome $workFasta"
    # IGV loads the genome asynchronously in batch mode; without a pause
    # the next command can run before the genome finishes loading.
    $lines += "setSleepInterval 2000"
    $lines += "new"
    $lines += "snapshotDirectory $workSnapDir"
    $lines += "maxPanelHeight 1000"

    foreach ($bam in $bamFiles) {
        $workBam = Join-Path $workBamDir $bam.Name
        Copy-Item -Path $bam.FullName -Destination $workBam -Force
        $baiSource = "$($bam.FullName).bai"
        if (Test-Path $baiSource) {
            Copy-Item -Path $baiSource -Destination "$workBam.bai" -Force
        }

        $baseName = [System.IO.Path]::GetFileNameWithoutExtension($bam.Name)
        $lines += "load $workBam"
        $lines += "goto $SequenceId"
        $lines += "snapshot ${baseName}_snapshot.png"
        $lines += "remove $workBam"
    }
    $lines += "exit"

    $bamCount = $bamFiles.Count

    if ($bamCount -gt 0) {
        Set-Content -Path $BatchFile -Value $lines -Encoding ASCII

        Write-Host "   Valid Reference Identified: $($faFile.Name)"
        Write-Host "   Sequence Targeted: [$SequenceId]"
        Write-Host "   Generating $bamCount snapshots using IGV..."

        Push-Location $IgvDir
        & "$JavaExec" -Xmx4000m `
            -Djava.net.useSystemProxies=true `
            -cp "lib\*" `
            org.broad.igv.ui.Main --batch "$BatchFile"
        Pop-Location

        # 7. Copy the resulting snapshots back to the real per-sample folder.
        $snaps = Get-ChildItem -Path $workSnapDir -Filter "*.png" -ErrorAction SilentlyContinue
        if ($snaps.Count -gt 0) {
            Copy-Item -Path $snaps.FullName -Destination $outDir -Force
            Write-Host "   Copied $($snaps.Count) snapshot(s) to $outDir"
        } else {
            Write-Host "   No snapshot images were produced - check igv.log in $IgvDir" -ForegroundColor Yellow
        }

        Remove-Item -Path $BatchFile -Force -ErrorAction SilentlyContinue
    } else {
        Write-Host "   No .bam files found inside $bamFolder" -ForegroundColor Yellow
    }

    # Clean up this sample's scratch folder
    Remove-Item -Path $workDir -Recurse -Force -ErrorAction SilentlyContinue
}

Write-Host "All directories completed!" -ForegroundColor Green
