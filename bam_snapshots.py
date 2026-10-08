#!/usr/bin/env python3
"""
bam_snapshots.py - headless, cross-platform BAM snapshot generator.

Replaces the IGV-based shell/PowerShell scripts. Same folder layout:

    main_dir/
        sample1/
            reference.fa        (any .fa/.fasta, ignoring names containing "_clean")
            BAM/*.bam           (a .bai index is optional)
        sample2/ ...

For every BAM a PNG (<bam name>_snapshot.png) is written to sample/Snapshots/.
Each PNG has a coverage track (bases that differ from the reference in >=20%
of reads are coloured) and a packed read pileup with mismatches, insertions
and deletions marked.

Requires: python>=3.8, numpy, matplotlib, and either pysam (fast, Mac/Linux)
or bamnostic (pure Python, works everywhere).  pip install numpy matplotlib bamnostic
"""
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

try:
    import pysam as _bam
    BACKEND = "pysam"
except ImportError:
    try:
        import bamnostic as _bam
        BACKEND = "bamnostic"
    except ImportError:
        sys.exit("Need pysam or bamnostic:  pip install bamnostic")

BASE_COLOURS = {"A": "#2ca02c", "C": "#1f4fd8", "G": "#e08a00", "T": "#d62728", "N": "#888888"}
M, I, D, N, S, H, P, EQ, X = range(9)  # BAM CIGAR op codes
REF_CONSUMING = {M, D, N, EQ, X}


def read_first_contig(fasta):
    """Return (name, sequence-as-uppercase-bytes-array) of the first record."""
    name, parts = None, []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
            elif name is not None:
                parts.append(line.strip().upper())
    if name is None:
        raise ValueError("no FASTA header found")
    return name, np.frombuffer("".join(parts).encode(), dtype="S1")


def analyse_bam(bam_path, contig, ref, start, end):
    """Single pass: coverage, per-base allele counts, and per-read geometry."""
    L = len(ref)
    cov_diff = np.zeros(L + 1, dtype=np.int64)
    counts = {b: np.zeros(L, dtype=np.int32) for b in "ACGT"}
    reads = []  # (start, end, reverse, mismatches[(pos, base)], insertions[pos], deletions[(s, e)])

    mode = "rb"
    with _bam.AlignmentFile(str(bam_path), mode) as af:
        try:
            it = af.fetch(contig, start, end)
        except Exception:  # no index: scan the whole file
            it = (r for r in af if r.reference_name == contig)
        for r in it:
            if r.is_unmapped or r.is_secondary or not r.cigartuples:
                continue
            qseq = (r.query_sequence or "").upper()
            q = np.frombuffer(qseq.encode(), dtype="S1") if qseq else None
            rpos, qpos = r.reference_start, 0
            mism, ins, dels = [], [], []
            for op, ln in r.cigartuples:
                if op in (M, EQ, X):
                    a, b = max(rpos, 0), min(rpos + ln, L)
                    if b > a:
                        cov_diff[a] += 1
                        cov_diff[b] -= 1
                        if q is not None:
                            seg = q[qpos + (a - rpos): qpos + (b - rpos)]
                            for base in "ACGT":
                                counts[base][a:b] += seg == base.encode()
                            bad = np.nonzero((seg != ref[a:b]) & (seg != b"N"))[0]
                            mism.extend((a + int(i), seg[i].decode()) for i in bad)
                    rpos += ln
                    qpos += ln
                elif op == I:
                    ins.append(rpos)
                    qpos += ln
                elif op == S:
                    qpos += ln
                elif op in (D, N):
                    if op == D:
                        dels.append((rpos, rpos + ln))
                        cov_diff[rpos] += 1  # deletions still count as spanned
                        cov_diff[min(rpos + ln, L)] -= 1
                    rpos += ln
            reads.append((r.reference_start, rpos, bool(r.is_reverse), mism, ins, dels))
    return np.cumsum(cov_diff[:L]), counts, reads


def pack_rows(reads, max_rows):
    """Greedy interval packing; reads that don't fit in max_rows are dropped."""
    row_ends, placed = [], []
    for rd in sorted(reads, key=lambda x: x[0]):
        for i, e in enumerate(row_ends):
            if rd[0] > e + 1:
                row_ends[i] = rd[1]
                placed.append((i, rd))
                break
        else:
            if len(row_ends) < max_rows:
                row_ends.append(rd[1])
                placed.append((len(row_ends) - 1, rd))
    return placed, len(row_ends)


def render(bam_path, out_png, contig, ref, region, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    start, end = region
    cov, counts, reads = analyse_bam(bam_path, contig, ref, start, end)
    n_total = len(reads)
    if n_total > args.max_reads:
        rng = np.random.default_rng(0)
        reads = [reads[i] for i in rng.choice(n_total, args.max_reads, replace=False)]
    placed, n_rows = pack_rows(reads, args.max_rows)
    n_rows = max(n_rows, 1)

    fig_h = 2.2 + n_rows * 0.09
    fig, (axc, axr) = plt.subplots(
        2, 1, figsize=(args.width, fig_h), sharex=True,
        gridspec_kw={"height_ratios": [2.0, max(n_rows * 0.09, 0.5)], "hspace": 0.04})

    xs = np.arange(len(ref))
    axc.fill_between(xs, cov, step="mid", color="#b8b8b8", linewidth=0)
    total = np.maximum(cov, 1)
    for base, colour in BASE_COLOURS.items():
        if base == "N":
            continue
        frac = counts[base] / total
        mism = (frac >= 0.2) & (ref != base.encode()) & (cov > 0)
        if mism.any():
            axc.vlines(xs[mism], 0, cov[mism], colors=colour, linewidth=1.2)
    axc.set_xlim(start, end)
    axc.set_ylim(0, max(int(cov[start:end].max()) if end > start else 1, 1) * 1.05)
    axc.set_ylabel("Coverage")
    axc.set_title(f"{Path(bam_path).name}   {contig}:{start + 1:,}-{end:,}   "
                  f"{n_total:,} reads" + (f" ({len(placed):,} shown)" if len(placed) < n_total else ""),
                  loc="left", fontsize=10)

    fwd, rev, ticks, tcol, insx, dels = [], [], [], [], [], []
    for row, (s, e, is_rev, mism, ins, dl) in placed:
        y = row
        poly = [(s, y + .1), (e, y + .1), (e, y + .9), (s, y + .9)]
        (rev if is_rev else fwd).append(poly)
        for p, b in mism:
            ticks.append([(p, y + .1), (p + 1, y + .1), (p + 1, y + .9), (p, y + .9)])
            tcol.append(BASE_COLOURS.get(b, "#888888"))
        insx.extend((p, y) for p in ins)
        dels.extend((a, b, y) for a, b in dl)
    axr.add_collection(PolyCollection(fwd, facecolor="#cfd8e8", edgecolor="none"))
    axr.add_collection(PolyCollection(rev, facecolor="#e8d0d0", edgecolor="none"))
    for a, b, y in dels:  # deletion = black line through the read
        axr.plot([a, b], [y + .5, y + .5], color="black", linewidth=1)
    if ticks:
        axr.add_collection(PolyCollection(ticks, facecolor=tcol, edgecolor=tcol, linewidth=0.4))
    if insx:
        axr.scatter(*zip(*[(p, y + .5) for p, y in insx]), marker="|", s=40, color="#7a2fbf", linewidths=1.2)
    axr.set_ylim(n_rows, 0)
    axr.set_yticks([])
    axr.set_xlabel(contig)
    axr.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    for ax in (axc, axr):
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)


def _job(task):
    bam, out, contig, ref, region, args = task
    try:
        render(bam, out, contig, ref, region, args)
        return bam, None
    except Exception as exc:  # report and keep going
        return bam, exc


def find_fasta(folder):
    for p in sorted(folder.iterdir()):
        if p.is_file() and p.suffix.lower() in (".fa", ".fasta") and "_clean" not in p.name.lower():
            return p
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("main_dir", type=Path)
    ap.add_argument("--region", help="START-END (1-based, on the first contig); default: whole contig")
    ap.add_argument("--width", type=float, default=14, help="figure width in inches")
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument("--max-reads", type=int, default=2000, help="reads drawn (coverage always uses all)")
    ap.add_argument("--max-rows", type=int, default=80, help="max stacked read rows")
    ap.add_argument("--jobs", type=int, default=None, help="parallel processes (default: CPU count)")
    args = ap.parse_args()

    if not args.main_dir.is_dir():
        sys.exit(f"Error: '{args.main_dir}' is not a directory.")
    print(f"Using {BACKEND} backend; processing {args.main_dir.resolve()}")

    tasks = []
    for folder in sorted(p for p in args.main_dir.iterdir() if p.is_dir()):
        fa, bam_dir = find_fasta(folder), folder / "BAM"
        if fa is None or not bam_dir.is_dir():
            print(f"  skip {folder.name}: needs a non-_clean .fa/.fasta and a BAM/ folder")
            continue
        try:
            contig, ref = read_first_contig(fa)
        except ValueError as exc:
            print(f"  skip {folder.name}: {exc}")
            continue
        region = (0, len(ref))
        if args.region:
            s, e = args.region.replace(",", "").split("-")
            region = (max(int(s) - 1, 0), min(int(e), len(ref)))
        out_dir = folder / "Snapshots"
        out_dir.mkdir(exist_ok=True)
        bams = sorted(bam_dir.glob("*.bam"))
        if not bams:
            print(f"  skip {folder.name}: no .bam files")
        for b in bams:
            tasks.append((b, out_dir / f"{b.stem}_snapshot.png", contig, ref, region, args))

    failures = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for bam, err in ex.map(_job, tasks):
            if err:
                failures += 1
                print(f"  FAILED {bam}: {err}")
            else:
                print(f"  ok {bam.parent.parent.name}/{bam.name}")
    print(f"Done: {len(tasks) - failures}/{len(tasks)} snapshots written.")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
