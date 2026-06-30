#!/usr/bin/env python
# coding: utf-8
"""
R4 track panels (a, d, g, h) for the novel-transcription figure.
Panel letters per the reordered MANUSCRIPT.md.

  a  7/2 x 7/4  schematic of the 3 antisense cases (spurious / read-through /
                embedded), genes in CASE_COLOR (b/c palette); antisense genes
                transparent + dashed; one gray isoform-span arrow at the bottom.
  d  7/2 x 7/4  his3 / MMSYN1_0918 antisense transcription (yeast HIS3 marker):
                gene track | + strand antisense isoform track | + strand depth.
  g  7/2 x 7/4  truly intergenic transcript between lap/0154 and 0155:
                gene track (both flanking genes) | isoform track | + strand depth.
  h  7/2 x 7/4  novel ORF NOVEL_PEP_002 (118 aa) in the 0591/0592 gap:
                gene/ORF track | isoform track (no depth, no AF3 inset).

Reads the canonical isoform table + GFF + PacBio plus-strand depth bedGraph.
Output: R4_panels/panel_{a,d,g,h}.pdf
Run from Syn1_Novel_ORF/.
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, FancyBboxPatch

mpl.rcParams.update({
    'font.size': 7, 'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Nimbus Sans', 'Helvetica', 'DejaVu Sans'],
    'axes.titlesize': 7, 'axes.labelsize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

HALF, QUART = 7 / 2, 7 / 4
OUT = 'R4_panels'
CHROM = 'CP002027.1'
ISO_TSV = '../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv'
GFF = '../Genomes_Input/syn1.genes.gff3'
DEPTH_PLUS = '../Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ.plus.bedGraph'
DEPTH_MINUS = '../Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ.minus.bedGraph'
GENOME_LEN = 1_078_809   # Syn1 (CP002027.1)
SYN3A_DEPTH_PLUS = '../Syn3A_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph/syn3A_rep1.plus.bedGraph'
SYN3A_DEPTH_MINUS = '../Syn3A_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph/syn3A_rep1.minus.bedGraph'
SYN3A_CHROM = 'CP016816.2'
SYN3A_LEN = 543_379

CASE_COLOR = {'spurious_prom': '#0072B2', 'read_through': '#D55E00', 'embedded': '#009E73'}

# ====================================================================== loaders
def load_genes():
    rows = []
    for line in open(GFF):
        if not line.strip() or line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) != 9 or p[2] != 'gene':
            continue
        a = dict(kv.split('=', 1) for kv in p[8].split(';') if '=' in kv)
        rows.append({'chrom': p[0], 'start0': int(p[3]) - 1, 'end0': int(p[4]), 'strand': p[6],
                     'locus_tag': a.get('locus_tag', ''),
                     'gene_name': a.get('Name') or a.get('gene') or a.get('locus_tag', ''),
                     'rna_type': a.get('rna_type', '')})
    return pd.DataFrame(rows)

GENES = load_genes()
ISO = pd.read_csv(ISO_TSV, sep='\t')

DEL_BED = '../Genome_Reduction/aln/raw/syn1_deleted_regions.bed'
def load_deletions():
    out = []
    for line in open(DEL_BED):
        if line.startswith('#') or not line.strip():
            continue
        p = line.split('\t')
        out.append((int(p[1]), int(p[2])))
    return out
DELETIONS = load_deletions()


def load_depth_window(win_s, win_e):
    """Per-base plus-strand depth over [win_s, win_e) via awk slice (fast)."""
    import subprocess
    cov = np.zeros(win_e - win_s)
    cmd = ['awk', '-F', '\t',
           f'$1=="{CHROM}" && $3>{win_s} && $2<{win_e}', DEPTH_PLUS]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    for line in out.splitlines():
        _, s, e, v = line.split('\t')
        s, e, v = int(s), int(e), float(v)
        a, b = max(s, win_s) - win_s, min(e, win_e) - win_s
        cov[a:b] = v
    return cov


_MEAN_TOTAL = [None]
def genome_mean_depth_total():
    """Genome-wide mean per-base (plus+minus) Syn1 PacBio depth (for x-mean normalisation)."""
    if _MEAN_TOTAL[0] is None:
        import subprocess
        tot = 0.0
        for path in (DEPTH_PLUS, DEPTH_MINUS):
            cmd = ['awk', '-F', '\t', f'$1=="{CHROM}"{{s+=($3-$2)*$4}} END{{print s+0}}', path]
            out = subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
            tot += float(out) if out else 0.0
        _MEAN_TOTAL[0] = tot / GENOME_LEN
    return _MEAN_TOTAL[0]


def _nice_top(m):
    """Smallest 'nice' tick >= m (1/1.5/2/3/5/7 x 10^k)."""
    if m <= 0:
        return 1.0
    mag = 10.0 ** np.floor(np.log10(m))
    for k in (1, 1.5, 2, 3, 5, 7, 10):
        if k * mag >= m - 1e-9:
            return k * mag
    return 10 * mag


def load_syn3a_depth_plus(q_s, q_e):
    """Per-base Syn3A Illumina plus-strand depth over [q_s, q_e) via awk slice."""
    import subprocess
    cov = np.zeros(q_e - q_s)
    cmd = ['awk', '-F', '\t', f'$1=="{SYN3A_CHROM}" && $3>{q_s} && $2<{q_e}', SYN3A_DEPTH_PLUS]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    for line in out.splitlines():
        _, s, e, v = line.split('\t')
        s, e, v = int(s), int(e), float(v)
        a, b = max(s, q_s) - q_s, min(e, q_e) - q_s
        cov[a:b] = v
    return cov


_SYN3A_MEAN_TOTAL = [None]
def syn3a_mean_depth_total():
    """Genome-wide mean per-base (plus+minus) Syn3A Illumina depth (for x-mean normalisation)."""
    if _SYN3A_MEAN_TOTAL[0] is None:
        import subprocess
        tot = 0.0
        for path in (SYN3A_DEPTH_PLUS, SYN3A_DEPTH_MINUS):
            cmd = ['awk', '-F', '\t', f'$1=="{SYN3A_CHROM}"{{s+=($3-$2)*$4}} END{{print s+0}}', path]
            out = subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
            tot += float(out) if out else 0.0
        _SYN3A_MEAN_TOTAL[0] = tot / SYN3A_LEN
    return _SYN3A_MEAN_TOTAL[0]


def gene_label(r):
    """gene_name/locusNum; if no gene_name, locusNum; pseudogenes without a name
    -> pseudo/locusNum (e.g. his3/0918, lap/0154, mmyCImod/0591, pseudo/0155, 0592)."""
    num = str(r.locus_tag).split('_')[-1]
    nm = str(r.gene_name).strip()
    if nm and nm != str(r.locus_tag) and nm.lower() != 'nan':
        return f'{nm}/{num}'
    if str(r.rna_type) == 'pseudo':
        return f'pseudo/{num}'
    return num


# ====================================================================== tx helpers
def tx(pos0, win_s, win_e, strand):
    return pos0 - win_s if strand == '+' else win_e - pos0

def tx_to_genome(x, win_s, win_e, strand):
    return win_s + x if strand == '+' else win_e - x

def _pack_rows(intervals, gap=20):
    order = sorted(range(len(intervals)), key=lambda i: intervals[i][0])
    ends, rows = [], [0] * len(intervals)
    for i in order:
        a, b = intervals[i]
        placed = False
        for r, e in enumerate(ends):
            if a >= e + gap:
                ends[r] = b; rows[i] = r; placed = True; break
        if not placed:
            ends.append(b); rows[i] = len(ends) - 1
    return rows


# ====================================================================== track drawers
def draw_gene_track(ax, win_s, win_e, strand, orf=None, label_genes=True):
    """Gene polygon-arrows (gray, pointing by gene strand) on a baseline; optional
    novel-ORF box (red). x is tx-space."""
    win_len = win_e - win_s
    ax.set_xlim(0, win_len); ax.set_ylim(0, 1.6)
    # syn3A deletion overlay (shaded band behind the genes)
    drew_del = False
    for d0, d1 in DELETIONS:
        if d1 <= win_s or d0 >= win_e:
            continue
        a, b = max(d0, win_s), min(d1, win_e)
        xa, xb = tx(a, win_s, win_e, strand), tx(b, win_s, win_e, strand)
        ax.axvspan(min(xa, xb), max(xa, xb), facecolor='#e8736a', alpha=0.17, lw=0, zorder=0)
        drew_del = True
    if drew_del:
        ax.text(0.01, 0.97, 'syn3A deletion (shaded)', transform=ax.transAxes,
                ha='left', va='top', fontsize=4.5, color='#c0392b')
    ax.hlines(0.6, 0, win_len, color='black', lw=1.0, zorder=1)
    g = GENES[(GENES.chrom == CHROM) & (GENES.end0 > win_s) & (GENES.start0 < win_e)]
    H, TRI = 0.30, 0.30
    for _, r in g.iterrows():
        x0, x1 = tx(int(r.start0), win_s, win_e, strand), tx(int(r.end0), win_s, win_e, strand)
        xl, xr = min(x0, x1), max(x0, x1)
        head = min(max(40, (xr - xl) * 0.2), xr - xl)
        right = (r.strand == strand)
        col = '#b0a0c8' if r.rna_type == 'pseudo' else '#7a7a7a'   # pseudogenes purple
        if right:
            tip, base = xr, xr - head
            v = [(xl, 0.6 - H/2), (base, 0.6 - H/2), (base, 0.6 - TRI/2), (tip, 0.6),
                 (base, 0.6 + TRI/2), (base, 0.6 + H/2), (xl, 0.6 + H/2)]
        else:
            tip, base = xl, xl + head
            v = [(xr, 0.6 - H/2), (base, 0.6 - H/2), (base, 0.6 - TRI/2), (tip, 0.6),
                 (base, 0.6 + TRI/2), (base, 0.6 + H/2), (xr, 0.6 + H/2)]
        anti = (r.strand != strand)                          # antisense to the shown transcription
        ax.add_patch(Polygon(v, closed=True, facecolor=col, edgecolor='black',
                             lw=0.5 if anti else 0.3, alpha=0.45 if anti else 1.0,
                             linestyle='--' if anti else '-', zorder=2))
        vis_l, vis_r = max(0, xl), min(win_len, xr)        # clipped span for label
        if label_genes and (vis_r - vis_l) > 0.04 * win_len:
            ax.text((vis_l + vis_r) / 2, 0.6 + H/2 + 0.08, gene_label(r), ha='center', va='bottom',
                    fontsize=5, color='#333', clip_on=True)
    if orf is not None:
        o0, o1, oname = orf
        x0, x1 = tx(o0, win_s, win_e, strand), tx(o1, win_s, win_e, strand)
        oxl, oxr = min(x0, x1), max(x0, x1)
        ax.add_patch(FancyBboxPatch((oxl, 1.15 - 0.16), oxr - oxl, 0.32,
                     boxstyle='round,pad=0', lw=0.8, edgecolor='#800000',
                     facecolor='#D6604D', alpha=0.9, zorder=3))
        ax.text((oxl + oxr) / 2, 1.15 + 0.20, oname, ha='center', va='bottom',
                fontsize=5, fontweight='bold', color='#800000', clip_on=True)
        for xv in (oxl, oxr):
            ax.axvline(xv, color='#800000', lw=0.6, ls='--', alpha=0.5, zorder=0)
    ax.set_yticks([]); ax.set_xticks([])
    ax.spines[['top', 'right', 'left', 'bottom']].set_visible(False)


def draw_isoform_track(ax, iso, win_s, win_e, strand, color='#1b6ca8', max_iso=70):
    win_len = win_e - win_s
    iso = iso.sort_values('n_reads', ascending=False).head(max_iso).copy()
    ints = []
    for _, r in iso.iterrows():
        a = tx(min(max(int(r.start0), win_s), win_e), win_s, win_e, strand)
        b = tx(min(max(int(r.end0), win_s), win_e), win_s, win_e, strand)
        ints.append((min(a, b), max(a, b)))
    rows = _pack_rows(ints)
    nmax = max(1, iso['n_reads'].max())
    for (xl, xr), ri, (_, r) in zip(ints, rows, iso.iterrows()):
        lw = float(np.clip(0.3 + 0.7 * np.log10(max(1, r['n_reads'])), 0.5, 3.0))
        al = 0.45 + 0.5 * (r['n_reads'] / nmax)
        ax.add_patch(FancyArrowPatch((xl, ri), (xr, ri), arrowstyle='-|>',
                     lw=lw, color=color, alpha=min(1.0, al), shrinkA=0, shrinkB=0,
                     mutation_scale=5, zorder=2))
    ax.set_xlim(0, win_len); ax.set_ylim(-1, (max(rows) if rows else 0) + 1.5)
    ax.set_yticks([])
    ax.set_ylabel('isoforms', fontsize=5)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.set_xticks([])


def draw_depth_track(ax, win_s, win_e, strand, color='#9ecae1'):
    cov = load_depth_window(win_s, win_e) / genome_mean_depth_total()   # x genome-mean (both strands)
    win_len = win_e - win_s
    xg = np.arange(win_len)
    if strand == '-':
        cov = cov[::-1]
    ax.fill_between(xg, 0, cov, color=color, lw=0, zorder=1)
    ax.plot(xg, cov, color='#3182bd', lw=0.4, zorder=2)
    m = float(cov.max())
    if m > 0:
        T = _nice_top(m)                                   # top tick = nice x-mean value
        ax.set_yticks([0, T])
        ax.set_yticklabels(['0', f'{T:.0f}×' if T >= 1 else f'{T:g}×'])
        ax.set_ylim(0, T * 1.02)
    else:
        ax.set_ylim(0, 1)
    ax.set_xlim(0, win_len)
    ax.set_ylabel('Depth\n(× mean)', fontsize=5)
    ax.tick_params(labelsize=5)
    ax.spines[['top', 'right']].set_visible(False)
    # genomic x ticks
    ticks = np.linspace(0, win_len, 5)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{tx_to_genome(int(t), win_s, win_e, strand)/1000:.1f}' for t in ticks])
    ax.set_xlabel('Syn1 Genome Position (kb)', fontsize=6)


def locus_panel(fname, win_s, win_e, strand, iso_sel, iso_color, orf=None, depth=True, seq=None):
    """Stacked gene | isoform | (depth) panel, born at 7/2 x 7/4. Optional protein
    `seq` is wrapped and placed in the blank right of the isoform track."""
    nrow = 3 if depth else 2
    hr = [1.0, 2.2, 1.1] if depth else [1.3, 2.4]
    fig, axes = plt.subplots(nrow, 1, figsize=(HALF, QUART), height_ratios=hr,
                             constrained_layout=True)
    draw_gene_track(axes[0], win_s, win_e, strand, orf=orf)
    draw_isoform_track(axes[1], iso_sel, win_s, win_e, strand, color=iso_color)
    if depth:
        draw_depth_track(axes[2], win_s, win_e, strand)
    else:
        draw_isoform_xaxis(axes[1], win_s, win_e, strand)
    if seq is not None:
        import textwrap
        axes[1].text(0.995, 0.55, textwrap.fill(seq, 24), transform=axes[1].transAxes,
                     ha='right', va='center', fontsize=5, family='monospace', color='#444')
    fig.savefig(f'{OUT}/{fname}', dpi=300); plt.close(fig)


def draw_isoform_xaxis(ax, win_s, win_e, strand):
    win_len = win_e - win_s
    ticks = np.linspace(0, win_len, 5)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{tx_to_genome(int(t), win_s, win_e, strand)/1000:.1f}' for t in ticks],
                       fontsize=5)
    ax.tick_params(labelsize=5)
    ax.spines['bottom'].set_visible(True)
    ax.set_xlabel('Syn1 Genome Position (kb)', fontsize=6)


# ====================================================================== panel a
def schematic_gene(ax, xl, xr, y, points_right, color, antisense, h=0.32):
    head = (xr - xl) * 0.22
    if points_right:
        tip, base = xr, xr - head
        v = [(xl, y - h/2), (base, y - h/2), (base, y - h/2 - 0.05), (tip, y),
             (base, y + h/2 + 0.05), (base, y + h/2), (xl, y + h/2)]
    else:
        tip, base = xl, xl + head
        v = [(xr, y - h/2), (base, y - h/2), (base, y - h/2 - 0.05), (tip, y),
             (base, y + h/2 + 0.05), (base, y + h/2), (xr, y + h/2)]
    ax.add_patch(Polygon(v, closed=True, facecolor=color,
                         alpha=0.35 if antisense else 0.92,
                         edgecolor=color, lw=1.0,
                         linestyle='--' if antisense else '-', zorder=2))


def panel_a():
    fig, ax = plt.subplots(figsize=(HALF, QUART), constrained_layout=True)
    cases = [('spurious_prom', 'Spurious promoter', 59, 2.0),
             ('read_through', 'Read-through', 26, 1.0),
             ('embedded', 'Embedded', 4, 0.0)]
    X0, X1 = 0.0, 10.0
    for key, label, n, y in cases:
        col = CASE_COLOR[key]
        ax.hlines(y, X0, X1, color='black', lw=0.9, zorder=1)          # genome line
        if key == 'spurious_prom':                                     # anti only
            schematic_gene(ax, 4.0, 6.4, y, points_right=False, color=col, antisense=True)
        elif key == 'read_through':                                    # sense -> anti
            schematic_gene(ax, 1.0, 4.0, y, points_right=True,  color=col, antisense=False)
            schematic_gene(ax, 5.6, 8.6, y, points_right=False, color=col, antisense=True)
        else:                                                          # sense -> anti -> sense
            schematic_gene(ax, 0.8, 3.0, y, points_right=True,  color=col, antisense=False)
            schematic_gene(ax, 3.6, 6.0, y, points_right=False, color=col, antisense=True)
            schematic_gene(ax, 6.6, 9.0, y, points_right=True,  color=col, antisense=False)
        ax.text(X0 - 0.3, y, f'{label}\n(n={n})', ha='right', va='center',
                fontsize=6, color=col)
    # isoform-span arrow drawn in the bottom margin (panel a's "x-axis"); the span
    # label is a real xlabel so its font matches panels b/c
    ax.annotate('', xy=(X1, -0.08), xytext=(X0, -0.08),
                xycoords=('data', 'axes fraction'), annotation_clip=False,
                arrowprops=dict(arrowstyle='-|>', color='#808080', lw=2.2,
                                shrinkA=0, shrinkB=0, mutation_scale=10))
    # cases stay at y=2/1/0; ylim -0.2..2.9 maps them to the b/c ridge-baseline
    # heights, so the three cases align row-for-row with panels b and c
    ax.set_xlim(-4.2, 10.3); ax.set_ylim(-0.2, 2.9)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel(r"Isoform span (5'$\to$3')", labelpad=18)
    fig.savefig(f'{OUT}/panel_a_schematic.pdf', dpi=300); plt.close(fig)


# ====================================================================== -10 box at novel TSS
def quantify_novel_promoters(out_txt=f'{OUT}/novel_promoter_minus10.txt'):
    """Scan the -10 promoter box at the TSS of the two non-canonical transcription
    events (his3/0918 antisense spurious promoter; intergenic lap-0154/0155), using
    the SAME algorithm as the canonical operon promoters (Syn1_Operon/promoter_motif.py,
    the single source). TSS = the 5' end (pos5p0) of the most-supported (+)-strand
    isoform in each case. Reports both -10 boxes (6-mer TANAAT, 9-mer TNNTANAAT) +
    the -35 region, so they're directly comparable to the canonical promoters."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Syn1_Operon'))
    import promoter_motif as pm

    his_s, his_e = 27639, 28301
    inter_s, inter_e = 197500, 201700
    cases = [
        ("his3/0918 antisense (spurious promoter)",
         ISO[(ISO.strand == '+') & (ISO.start0 < his_e) & (ISO.end0 > his_s) & (ISO.n_reads >= 10)]),
        ("intergenic lap/0154 - 0155",
         ISO[(ISO.strand == '+') & (ISO.start0 < inter_e) & (ISO.end0 > inter_s) &
             (ISO.n_reads >= 10) & (ISO.frac_intergenic > 0.5)]),
    ]
    L = ["-10 PROMOTER BOX AT NON-CANONICAL TSS",
         "=" * 52,
         "TSS = 5' end (pos5p0) of the most-supported (+)-strand isoform; scan uses",
         "promoter_motif.scan_minus10 (identical algorithm to the canonical operons).",
         "Canonical baseline: TANAAT 6-mer 87/127 (69%); extended TNNTANAAT 52/127 (41%).",
         ""]
    for label, sel in cases:
        if len(sel) == 0:
            L += [f"{label}: no isoforms passing the filter", ""]; continue
        top = sel.sort_values('n_reads', ascending=False).iloc[0]
        tss0 = int(top['pos5p0'])
        r = pm.scan_minus10(tss0, CHROM, '+')
        m35 = pm.extract_tx_kmer(tss0, CHROM, '+', -37, -31)
        at35 = sum(c in 'AT' for c in m35)
        L += [
            label,
            f"  isoforms (n_reads>=10): {len(sel)};  top {top['isoform_id']} "
            f"n_reads={int(top['n_reads'])};  TSS(pos5p0)={tss0}",
            f"  -10 6-mer (TANAAT,   [-12,-7]): {str(r['minus10_6mer']):>9}  "
            f"shift {r['shift6']:+d}  match={r['match6']}  mm={r['mm6']}",
            f"  -10 9-mer (TNNTANAAT,[-15,-7]): {str(r['minus10_9mer']):>9}  "
            f"shift {r['shift9']:+d}  match={r['match9']}  mm={r['mm9']}",
            f"  -35 region          [-37,-31]: {str(m35):>9}  (AT {at35}/{len(m35) or 0})",
            f"  -10 tier: {r['motif_tier']}",
            "",
        ]
    txt = "\n".join(L)
    with open(out_txt, 'w') as fh:
        fh.write(txt + "\n")
    print(txt)
    print(f"Saved: {out_txt}")


def panel_d_his3():
    """his3/0918 antisense over-transcription, Syn1 vs Syn3A, on a deletion-junction
    relative axis: rel = syn1_pos - deletion_end (his3 positive, deleted upstream negative).
    Tracks: genes | Syn1 + isoforms | Syn1 + depth (x-mean, blue) | Syn3A + depth (x-mean, red)."""
    D1, D3 = 27638, 18715              # junction: Syn1 deletion end <-> Syn3A his3 start0
    win_s, win_e = 26900, 29000
    win_len = win_e - win_s
    his_s, his_e = 27639, 28301
    sel = ISO[(ISO.strand == '+') & (ISO.start0 < his_e) & (ISO.end0 > his_s) & (ISO.n_reads >= 10)]

    fig, axes = plt.subplots(4, 1, figsize=(HALF, QUART),
                             height_ratios=[1.0, 0.75, 1.0, 1.0], constrained_layout=True)
    draw_gene_track(axes[0], win_s, win_e, '+')
    draw_isoform_track(axes[1], sel, win_s, win_e, '+', color=CASE_COLOR['spurious_prom'], max_iso=10)
    draw_depth_track(axes[2], win_s, win_e, '+')                      # Syn1 PacBio + depth (blue)
    axes[2].set_xticks([]); axes[2].set_xlabel('')
    axes[2].set_ylabel('Syn1\n(× mean)', fontsize=5, color='#3182bd')

    # Syn3A Illumina + depth, aligned by the junction, normalised to its own mean (red)
    q_s = D3 + (win_s - D1)
    cov3 = load_syn3a_depth_plus(q_s, q_s + win_len) / syn3a_mean_depth_total()
    ax3, xg = axes[3], np.arange(win_len)
    ax3.fill_between(xg, 0, cov3, color='#f3b0ad', lw=0, zorder=1)
    ax3.plot(xg, cov3, color='#c0392b', lw=0.4, zorder=2)
    m = float(cov3.max())
    if m > 0:
        T = _nice_top(m)
        ax3.set_yticks([0, T]); ax3.set_yticklabels(['0', f'{T:.0f}×' if T >= 1 else f'{T:g}×'])
        ax3.set_ylim(0, T * 1.02)
    else:
        ax3.set_ylim(0, 1)
    ax3.set_xlim(0, win_len)
    ax3.set_ylabel('Syn3A\n(× mean)', fontsize=5, color='#c0392b')
    ax3.tick_params(labelsize=5)
    ax3.spines[['top', 'right']].set_visible(False)
    off = D1 - win_s                                                  # rel -> tx: tx = rel + off
    rel_ticks = [-500, 0, 500, 1000]
    ax3.set_xticks([r + off for r in rel_ticks])
    ax3.set_xticklabels([str(r) for r in rel_ticks], fontsize=5)
    ax3.set_xlabel('Relative genome position (bp)', fontsize=6)

    fig.savefig(f'{OUT}/panel_d_his3_antisense.pdf', dpi=300); plt.close(fig)


# ====================================================================== main
def main():
    panel_a()

    # d: his3 / 0918 antisense over-transcription -- Syn1 vs Syn3A, deletion-junction axis
    panel_d_his3()

    # g: intergenic transcript between lap/0154 and 0155 (+ strand)
    win_s, win_e = 197500, 201700
    sel = ISO[(ISO.strand == '+') & (ISO.start0 < win_e) & (ISO.end0 > win_s) &
              (ISO.n_reads >= 10) & (ISO.frac_intergenic > 0.5)]
    locus_panel('panel_g_intergenic_0154_0155.pdf', win_s, win_e, '+', sel,
                iso_color='#1b9e77', depth=True)

    # h: novel ORF NOVEL_PEP_002 (- strand), gene + isoform only, with the 118-aa sequence
    orf_s, orf_e = 728399, 728756
    win_s, win_e = 726400, 730800
    PEP002 = ("MNTKQDYINKIDNLVLNSELNYDQKNLIISILNKFDDNDINLHNVYQFLIKRVKLGFTFDIAPSVDSDQVAI"
              "LSKDETRSFNNNKNNNNILIIGENYDALKNLIVAERERERERAGRC")
    sel = ISO[(ISO.strand == '-') & (ISO.start0 <= orf_s) & (ISO.end0 >= orf_e) & (ISO.n_reads >= 10)]
    locus_panel('panel_h_novel_orf_0592.pdf', win_s, win_e, '-', sel,
                iso_color='#1b6ca8', orf=(orf_s, orf_e, 'NOVEL_PEP_002 (118 aa)'),
                depth=False, seq=PEP002)

    print(f"Saved panels a, d, g, h to {OUT}/")
    quantify_novel_promoters()


if __name__ == '__main__':
    main()
