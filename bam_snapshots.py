#!/usr/bin/env python3
"""
bam_snapshots.py - headless, cross-platform BAM snapshot generator.

Expected folder layout:

    main_dir/
        sample1/
            reference.fa        (any .fa/.fasta, ignoring names containing "_clean")
            BAM/*.bam           (a .bai index is optional)
            *.gbk               (optional: annotated map -> "features" track)
            **/*.vcf            (optional: variant calls -> "variants" track;
                                 if absent, variants are called from the BAM)
        sample2/ ...

For every BAM, in sample/Snapshots/:

    <bam>_overview.png   whole plasmid: coverage, variants, features, read stack
    <bam>_zoom_<n>.png   base-level views around the highest-frequency variants
                         (or the lowest-coverage spot if there are none)

Runs on Windows, Mac and Linux. BAMs are processed in parallel.
Requires: python>=3.8, numpy, Pillow>=10.1 and either pysam (fast, Mac/Linux)
or bamnostic (pure Python, works everywhere):  pip install -r requirements.txt
"""
import argparse
import random
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

try:
    import pysam as _bam
    BACKEND = "pysam"
except ImportError:
    try:
        import bamnostic as _bam
        BACKEND = "bamnostic"
    except ImportError:
        sys.exit("Need pysam or bamnostic:  pip install bamnostic")

BASE_COL = {"A": (46, 160, 67), "C": (31, 111, 235), "G": (219, 134, 18), "T": (207, 34, 46), "N": (150, 150, 150)}
TITLE_COL = (0, 85, 150)
GREY = (120, 128, 138)
LIGHT = (232, 236, 241)
FEATURE_COL = {"CDS": (112, 144, 196), "promoter": (115, 176, 120), "terminator": (200, 120, 120),
               "rep_origin": (224, 170, 80), "primer_bind": (160, 160, 160), "misc_feature": (170, 170, 190)}
M, I, D, N, S, H, P, EQ, X = range(9)  # BAM CIGAR op codes
SKIP_FLAGS = 4 | 256 | 1024 | 2048     # unmapped, secondary, duplicate, supplementary


def font(size):
    try:
        return ImageFont.load_default(size=size)
    except TypeError:  # Pillow < 10.1
        return ImageFont.load_default()


# ---------------------------------------------------------------- inputs

def read_ref(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
            elif name is not None:
                seq.append(line.strip().upper())
    return name or "ref", "".join(seq)


def load_variants_vcf(folder):
    """VarScan-style VCFs under folder -> [{pos, ref, alt, freq(0-100), dp, kind}], deduplicated."""
    vcfs = [p for p in folder.rglob("*.vcf") if "Snapshots" not in p.parts]
    merged = [p for p in vcfs if "merged" in p.name.lower()]
    seen, out = set(), []
    for vcf in sorted(merged or vcfs):
        for line in vcf.read_text(errors="ignore").splitlines():
            f = line.split("\t")
            if line.startswith("#") or len(f) < 10 or (f[1], f[3], f[4]) in seen:
                continue
            seen.add((f[1], f[3], f[4]))
            d = dict(zip(f[8].split(":"), f[9].split(":")))
            try:
                freq = float(str(d.get("FREQ", "0")).strip("%"))
            except ValueError:
                freq = 0.0
            try:
                dp = int(d.get("DP", d.get("SDP", 0)) or 0)
            except ValueError:
                dp = 0
            out.append({"pos": int(f[1]), "ref": f[3], "alt": f[4], "freq": freq, "dp": dp,
                        "kind": "indel" if len(f[3]) != len(f[4]) else "snp"})
    return sorted(out, key=lambda v: v["pos"])


def load_features(folder):
    """Minimal GenBank feature-table parser -> [{type, start, end, label}] (1-based, inclusive)."""
    gbks = list(folder.glob("*.gbk")) + list(folder.glob("*.gb")) + list(folder.glob("*.gbff"))
    if not gbks:
        return []
    gbk = sorted(gbks, key=lambda p: "annotated" not in p.name.lower())[0]
    feats, cur = [], None
    for line in gbk.read_text(errors="ignore").splitlines():
        m = re.match(r"^ {5}(\S+)\s+(?:complement\()?(?:join\()?<?(\d+)\.\.>?(\d+)", line)
        if m:
            cur = {"type": m.group(1), "start": int(m.group(2)), "end": int(m.group(3)), "label": ""}
            if cur["type"] not in ("source", "gene"):
                feats.append(cur)
            continue
        m = re.match(r'^\s+/(label|gene|product)="?([^"]*)', line)
        if m and cur is not None and not cur["label"]:
            cur["label"] = m.group(2)
    return feats


# ---------------------------------------------------------------- BAM parsing

class Read:
    __slots__ = ("start", "end", "mism", "dels", "ins", "rev", "blocks")

    def __init__(self, start, end, mism, dels, ins, rev, blocks):
        self.start, self.end, self.mism, self.dels, self.ins, self.rev, self.blocks = \
            start, end, mism, dels, ins, rev, blocks


def parse_read(r, ref_arr, min_bq):
    """pysam/bamnostic alignment -> Read (None if unusable). Mismatches below min_bq are ignored."""
    if r.flag & SKIP_FLAGS or not r.cigartuples:
        return None
    L = len(ref_arr)
    q = np.frombuffer((r.query_sequence or "").upper().encode(), dtype="S1")
    qual = r.query_qualities
    qual = np.asarray(qual) if qual is not None and len(qual) else None
    rpos, qpos = r.reference_start, 0
    mism, dels, ins, blocks = [], [], [], []
    for op, ln in r.cigartuples:
        if op in (M, EQ, X):
            a, b = max(rpos, 0), min(rpos + ln, L)
            if b > a:
                blocks.append((a, b))
                if len(q):
                    seg = q[qpos + (a - rpos): qpos + (b - rpos)]
                    bad = (seg != ref_arr[a:b]) & (seg != b"N")
                    if qual is not None:
                        bad &= qual[qpos + (a - rpos): qpos + (b - rpos)] >= min_bq
                    mism.extend((a + int(i), seg[i].decode()) for i in np.nonzero(bad)[0])
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
            rpos += ln
    return Read(r.reference_start, rpos, mism, dels, ins, bool(r.flag & 16), blocks)


def iter_alignments(af, contig, start, end):
    try:
        yield from af.fetch(contig, start, end)  # needs an index
    except Exception:
        for r in af:                             # no index: scan everything
            if r.reference_name == contig and r.reference_start < end and (r.reference_end or 0) > start:
                yield r


def sample_reads(bam_path, contig, ref_arr, start, end, min_bq, cap, seed=1):
    """Reservoir-sample up to cap reads overlapping [start, end). Returns (reads, n_total)."""
    rng, reads, n = random.Random(seed), [], 0
    with _bam.AlignmentFile(str(bam_path), "rb") as af:
        for r in iter_alignments(af, contig, start, end):
            rd = parse_read(r, ref_arr, min_bq)
            if rd is None:
                continue
            n += 1
            if len(reads) < cap:
                reads.append(rd)
            else:
                j = rng.randrange(n)
                if j < cap:
                    reads[j] = rd
    reads.sort(key=lambda x: x.start)
    return reads, n


def scan_bam(bam_path, contig, ref_arr, min_bq, cap, seed=1):
    """One pass over the whole contig: depth, alt-allele/indel counts, and a read sample."""
    L = len(ref_arr)
    diff = np.zeros(L + 1, dtype=np.int64)
    alt, ins_n, del_n = {}, {}, {}
    rng, reads, n = random.Random(seed), [], 0
    with _bam.AlignmentFile(str(bam_path), "rb") as af:
        for r in iter_alignments(af, contig, 0, L):
            rd = parse_read(r, ref_arr, min_bq)
            if rd is None:
                continue
            n += 1
            for a, b in rd.blocks:
                diff[a] += 1
                diff[b] -= 1
            for k in rd.mism:
                alt[k] = alt.get(k, 0) + 1
            for p in set(rd.ins):
                ins_n[p] = ins_n.get(p, 0) + 1
            for a, b in rd.dels:
                del_n[a] = del_n.get(a, 0) + 1
            if len(reads) < cap:
                reads.append(rd)
            else:
                j = rng.randrange(n)
                if j < cap:
                    reads[j] = rd
    reads.sort(key=lambda x: x.start)
    return np.cumsum(diff[:L]), alt, ins_n, del_n, reads, n


def call_variants(ref_seq, depth, alt, ins_n, del_n, min_freq=0.05, min_reads=2, min_dp=10):
    """Simple frequency-threshold caller used when no VCF is available."""
    out = []
    for (p, b), c in alt.items():
        dp = int(depth[p])
        if c >= min_reads and dp >= min_dp and c / dp >= min_freq:
            out.append({"pos": p + 1, "ref": ref_seq[p], "alt": b, "freq": 100.0 * c / dp, "dp": dp, "kind": "snp"})
    for tag, table in (("ins", ins_n), ("del", del_n)):
        for p, c in table.items():
            dp = int(depth[min(p, len(depth) - 1)])
            if c >= min_reads and dp >= min_dp and c / dp >= min_freq:
                out.append({"pos": p + 1, "ref": ref_seq[p], "alt": tag, "freq": 100.0 * c / dp, "dp": dp,
                            "kind": "indel"})
    return sorted(out, key=lambda v: v["pos"])


# ---------------------------------------------------------------- drawing

def pack_rows(reads, max_rows, gap=2):
    rows, ends = [], []
    for r in reads:
        for i, e in enumerate(ends):
            if r.start > e + gap:
                rows[i].append(r)
                ends[i] = r.end
                break
        else:
            if len(rows) < max_rows:
                rows.append([r])
                ends.append(r.end)
    return rows


def nice_step(span, target=8):
    raw = max(1, span / target)
    for s in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000):
        if s >= raw:
            return s
    return 100000


def draw_ruler(d, x0, x1, y, start, end, f):
    span, step = end - start, nice_step(end - start)
    d.line([(x0, y), (x1, y)], fill=GREY)
    for p in range((start // step + 1) * step, end, step):
        x = x0 + (p - start) / span * (x1 - x0)
        d.line([(x, y), (x, y + 5)], fill=GREY)
        d.text((x + 3, y + 4), f"{p / 1000:g} kb" if step >= 500 else f"{p}", fill=GREY, font=f)


def draw_features(d, features, x0, x1, y, h, start, end, f):
    span = end - start
    for ft in features:
        s, e = max(ft["start"] - 1, start), min(ft["end"], end)
        if e <= s:
            continue
        xa = x0 + (s - start) / span * (x1 - x0)
        xb = x0 + (e - start) / span * (x1 - x0)
        d.rectangle([xa, y, max(xb, xa + 1), y + h], fill=FEATURE_COL.get(ft["type"], (170, 170, 190)))
        if xb - xa > 46:
            d.text((xa + 3, y + 1), str(ft["label"])[:int((xb - xa) / 7)], fill=(255, 255, 255), font=f)


def draw_variants(d, variants, x0, x1, y, h, start, end):
    span = end - start
    for v in variants:
        if not (start < v["pos"] <= end):
            continue
        x = x0 + (v["pos"] - 1 - start) / span * (x1 - x0)
        col = (31, 111, 235) if v["kind"] == "snp" else (219, 134, 18)
        hh = h * (0.45 + 0.55 * min(1.0, v["freq"] / 100.0))
        d.rectangle([x - 1, y + h - hh, x + 1, y + h], fill=col)


def render(reads, n_reads, ref_seq, depth, features, variants, start, end, title, out,
           width=1500, rows_max=60, row_h=4, base_level=False, mismatch_positions=None):
    """Draw the window [start, end) (0-based, end exclusive)."""
    start, end = max(0, start), min(len(ref_seq), end)
    span = end - start
    f, fs, fb = font(12), font(11), font(15)
    rows = pack_rows(reads, rows_max)
    x0, x1 = 70, width - 20
    px_per_bp = (x1 - x0) / span
    cov_h, var_h, feat_h = 70, 22, 18
    seq_h = 18 if base_level else 0
    y_cov = 30
    y_var = y_cov + cov_h + 6
    y_feat = y_var + var_h + 4
    y_ruler = y_feat + feat_h + 8
    y_seq = y_ruler + 26
    y_reads = y_seq + seq_h + 6
    height = y_reads + max(len(rows), 1) * row_h + 34
    img = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(img)
    d.text((x0, 6), title, fill=TITLE_COL, font=fb)

    # coverage
    win = depth[start:end]
    dmax = max(int(win.max()), 1) if len(win) else 1
    d.rectangle([x0, y_cov, x1, y_cov + cov_h], fill=(247, 249, 251), outline=LIGHT)
    cols = max(1, int(x1 - x0))
    for cx in range(cols):
        a = start + int(cx / cols * span)
        b = max(a + 1, start + int((cx + 1) / cols * span))
        seg = depth[a:b]
        v = float(seg.mean()) if len(seg) else 0
        d.line([(x0 + cx, y_cov + cov_h), (x0 + cx, y_cov + cov_h - v / dmax * (cov_h - 2))], fill=(120, 160, 205))
    d.text((6, y_cov), f"{dmax}x", fill=GREY, font=fs)
    d.text((6, y_cov + cov_h - 12), "0", fill=GREY, font=fs)
    d.text((6, y_cov + cov_h // 2 - 6), "depth", fill=GREY, font=fs)

    d.text((6, y_var + 4), "variants", fill=GREY, font=fs)
    draw_variants(d, variants, x0, x1, y_var, var_h, start, end)
    d.text((6, y_feat + 2), "features", fill=GREY, font=fs)
    d.rectangle([x0, y_feat, x1, y_feat + feat_h], fill=(247, 249, 251))
    draw_features(d, features, x0, x1, y_feat, feat_h, start, end, fs)
    draw_ruler(d, x0, x1, y_ruler, start, end, fs)

    if base_level:
        for i in range(span):
            b = ref_seq[start + i]
            x = x0 + i * px_per_bp
            d.rectangle([x, y_seq, x + px_per_bp - 1, y_seq + seq_h - 2], fill=BASE_COL.get(b, GREY))
            d.text((x + px_per_bp / 2 - 4, y_seq), b, fill="white", font=f)

    d.text((6, y_reads), "reads", fill=GREY, font=fs)
    for ri, row in enumerate(rows):
        y = y_reads + ri * row_h
        for r in row:
            xa = x0 + (max(r.start, start) - start) * px_per_bp
            xb = x0 + (min(r.end, end) - start) * px_per_bp
            if xb <= xa:
                continue
            d.rectangle([xa, y, xb, y + row_h - 2], fill=(205, 212, 222))
            for ds, de in r.dels:
                xs = x0 + (max(ds, start) - start) * px_per_bp
                xe = x0 + (min(de, end) - start) * px_per_bp
                if xe > xs:
                    d.rectangle([xs, y + row_h // 2 - 1, xe, y + row_h // 2], fill=(20, 20, 20))
            for ip in r.ins:
                if start <= ip < end:
                    xi = x0 + (ip - start) * px_per_bp
                    d.rectangle([xi - 1, y - 1, xi + 1, y + row_h - 1], fill=(130, 40, 160))
            for p, b in r.mism:
                if (mismatch_positions is not None and p not in mismatch_positions) or not (start <= p < end):
                    continue
                xm = x0 + (p - start) * px_per_bp
                d.rectangle([xm, y, xm + max(px_per_bp, 2.0), y + row_h - 2], fill=BASE_COL.get(b, GREY))
                if base_level and px_per_bp >= 9 and row_h >= 8:
                    d.text((xm + 2, y - 1), b, fill="white", font=fs)
    shown = sum(len(r) for r in rows)
    d.text((x0, height - 22),
           f"{n_reads:,} reads in view (duplicates hidden); drawing {shown:,}.  Coloured marks = bases that differ "
           f"from the reference (A green, C blue, G orange, T red); black = deletion; purple = insertion.",
           fill=GREY, font=fs)
    img.save(out, optimize=True)


def pick_sites(variants, depth, max_zoom, window):
    sites = []
    for v in sorted(variants, key=lambda v: (-v["freq"], -v["dp"])):
        if v["freq"] < 20 or v["dp"] < 10:
            continue
        if all(abs(v["pos"] - s) > window for s, _ in sites):
            sites.append((v["pos"], f"{v['kind'].upper()} at {v['pos']}: {v['ref']}>{v['alt']} ({v['freq']:.0f}%)"))
        if len(sites) >= max_zoom:
            break
    if not sites and len(depth):
        low = int(np.argmin(depth))
        if depth[low] < 10:
            sites.append((low + 1, f"Lowest coverage point: position {low + 1} ({int(depth[low])}x)"))
    return sites


# ---------------------------------------------------------------- driver

def process_bam(task):
    bam, sample_dir, out_dir, contig, ref_seq, args = task
    try:
        ref_arr = np.frombuffer(ref_seq.encode(), dtype="S1")
        stem = Path(bam).stem
        depth, alt, ins_n, del_n, reads, n = scan_bam(bam, contig, ref_arr, args.min_bq, args.max_reads)
        features = load_features(sample_dir)
        variants = load_variants_vcf(sample_dir) or call_variants(
            ref_seq, depth, alt, ins_n, del_n, args.min_var_freq)
        # overview: only colour mismatches at variant positions, so random errors don't look like noise
        var_pos = {v["pos"] - 1 + k for v in variants for k in range(max(len(v["ref"]), len(v["alt"])))}
        render(reads, n, ref_seq, depth, features, variants, 0, len(ref_seq),
               f"{stem}  |  whole plasmid ({len(ref_seq):,} bp)", out_dir / f"{stem}_overview.png",
               rows_max=70, row_h=3, mismatch_positions=var_pos)
        for i, (pos, label) in enumerate(pick_sites(variants, depth, args.max_zoom, args.zoom_window), start=1):
            s = max(0, pos - 1 - args.zoom_window // 2)
            e = min(len(ref_seq), s + args.zoom_window)
            zr, zn = sample_reads(bam, contig, ref_arr, s, e, args.min_bq, 2500)
            render(zr, zn, ref_seq, depth, features, variants, s, e, f"{stem}  |  {label}",
                   out_dir / f"{stem}_zoom_{i}.png", rows_max=45, row_h=9, base_level=True)
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
    ap.add_argument("--max-zoom", type=int, default=3, help="zoom PNGs per BAM")
    ap.add_argument("--zoom-window", type=int, default=100, help="bp shown in each zoom")
    ap.add_argument("--max-reads", type=int, default=6000, help="reads drawn in the overview")
    ap.add_argument("--min-bq", type=int, default=20, help="ignore mismatches with base quality below this")
    ap.add_argument("--min-var-freq", type=float, default=0.05,
                    help="frequency for BAM-called variants (used only when no VCF is found)")
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
        contig, ref_seq = read_ref(fa)
        if not ref_seq:
            print(f"  skip {folder.name}: empty reference")
            continue
        out_dir = folder / "Snapshots"
        out_dir.mkdir(exist_ok=True)
        bams = sorted(bam_dir.glob("*.bam"))
        if not bams:
            print(f"  skip {folder.name}: no .bam files")
        tasks += [(b, folder, out_dir, contig, ref_seq, args) for b in bams]

    failures = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for bam, err in ex.map(process_bam, tasks):
            if err:
                failures += 1
                print(f"  FAILED {bam}: {err!r}")
            else:
                print(f"  ok {bam.parent.parent.name}/{bam.name}")
    print(f"Done: {len(tasks) - failures}/{len(tasks)} BAMs snapshotted.")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
