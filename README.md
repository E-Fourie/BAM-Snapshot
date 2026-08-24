# IGV Batch Snapshot Scripts

Two versions of the same tool — pick the one for your OS. Both scripts expect
the same layout: a parent folder containing one subfolder per sample, each
with a reference fasta (any `.fa`/`.fasta`, excluding anything with
`_clean` in the name) and a `BAM/` subfolder. Snapshots are written to a
`Snapshots/` subfolder created inside each sample folder.

## Mac or Linux — `igv_batch_snapshots.sh`

Place it next to your extracted IGV folder (the one containing `lib/`), then:

```bash
chmod +x igv_batch_snapshots.sh
./igv_batch_snapshots.sh /path/to/your/main/directory
```

Auto-detects Mac vs Linux and adjusts the `find` syntax and Java lookup
accordingly.

## Windows — `igv_batch_snapshots.ps1`

Place it next to your extracted IGV folder (the one containing `lib`), then
either:

- Right-click the file → **Run with PowerShell**, or
- From a PowerShell prompt:

```powershell
.\igv_batch_snapshots.ps1 "C:\path\to\your\main\directory"
```

If Windows blocks the script with an execution-policy error, run this once
in the same PowerShell window first, then re-run the command above:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

## Java

Both scripts look for a Java 21 runtime in this order:

1. A JDK bundled inside the IGV folder itself (what IGV's own downloads
   ship with) — preferred, avoids version mismatches.
2. Common install locations (Homebrew on Mac; Adoptium/Temurin/Oracle on
   Windows).
3. Whatever `java`/`java.exe` is on your PATH.

If none are found, install Java 21 (e.g. Eclipse Temurin) and try again.

## Notes

- **The actual cause of the `GenomeManager.getCurrentGenome() is null` error:**
  IGV's batch-file parser has a long-standing, unfixed bug
  ([igvteam/igv#1258](https://github.com/igvteam/igv/issues/1258)) — it does
  **not** strip quote characters from quoted paths. Quoting a path (which is
  the usual workaround for spaces) makes IGV search for a filename that
  literally includes the quote marks, fail to find it, and silently fail to
  load the genome. The very next command then hits a null genome and throws
  that exact error — 100% reproducibly, even on paths with no spaces at all,
  since quoting itself is what breaks it.

  There's no working quoting syntax for this in IGV, so both scripts now
  avoid the problem entirely: every fasta/BAM/index file IGV needs is copied
  into a throwaway, guaranteed-space-free scratch folder first
  (`$TMPDIR/igv_batch.XXXXXX` on Mac/Linux, `C:\ProgramData\PSBA_IGV_scratch\<guid>`
  on Windows), and the batch file references *those* paths, completely
  unquoted. Resulting snapshots are copied back into the real per-sample
  `Snapshots/` folder afterwards, and the scratch folder is deleted.

- **A second, distinct cause of the same `getCurrentGenome() is null` error:**
  the batch file's `genome` line must come *before* the `new` line, not
  after. IGV's `new` command calls `newSession()`, which resets the view
  using whatever genome is currently active — and on a machine with no
  cached "last used genome" (e.g. a fresh IGV preferences directory, or the
  first time IGV has ever been run there), nothing is loaded yet when the
  batch script starts. `new` then throws this exact NullPointerException
  before the script's own `genome` command ever gets a chance to run. Both
  scripts now write `genome` first, so this no longer depends on IGV having
  a cached genome from some earlier session.

- Both scripts set `setSleepInterval 2000` at the start of each batch file.
  IGV loads the genome **asynchronously** in batch mode — without a pause,
  the next command (`load`/`goto`) can fire before the genome finishes
  loading. If you still see occasional issues on a slow disk or large
  fasta, try bumping the 2000 (ms) to something higher, e.g. 4000.
