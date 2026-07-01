#!/usr/bin/env python
"""
R5_figure_panels.py  --  born-at-size panels for Results figure R5
("Genome minimization remodels operon structure").

Pattern: one function per panel, each READS the already-computed analysis tables
(chiefly the SI workbook genome_reduction.xlsx) and emits one PDF into R5_panels/
at final print size per OUTPUT.md (Arial 5-7 pt, pdf.fonttype 42, no title,
constrained_layout). Recomputes nothing; statistics mirror 09/10.

Panels (R5 caption):
  a  genome map: 95 deletions overlaid on the operon map, by deletion class   [TODO]
  b  fusion junction DEL_014 (OP_00043 -> OP_00050)                           [TODO]
  c  transcript fold change by gene_impact_class (promoter_lost robustly down) [this]
  d  HupA operon, decapitated (promoter inside deleted MMSYN1_0349)           [TODO]

Run:  python R5_figure_panels.py c        # one panel
      python R5_figure_panels.py          # all implemented panels
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.size": 7,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

GR        = os.path.dirname(os.path.abspath(__file__))
SI_XLSX   = os.path.join(GR, "genome_reduction.xlsx")
OUTDIR    = os.path.join(GR, "R5_panels")
STATS_TXT = os.path.join(OUTDIR, "R5_panel_stats.txt")
os.makedirs(OUTDIR, exist_ok=True)

_log = []
def log(msg=""):
    print(msg)
    _log.append(msg)

# gene_impact_class ordered by predicted transcriptional effect: down -> neutral -> up
CLASS_ORDER = ["promoter_lost", "promoter_disconnected", "promoter_proximity_changed",
               "context_only", "unaffected", "readthrough_exposed", "new_promoter_fusion"]
CLASS_LABEL = {
    "promoter_lost":              "Promoter\nlost",
    "promoter_disconnected":      "Promoter\ndisconn.",
    "promoter_proximity_changed": "Proximity\nchanged",
    "context_only":               "Context\nonly",
    "unaffected":                 "Unaffected",
    "readthrough_exposed":        "Readthrough\nexposed",
    "new_promoter_fusion":        "New\nfusion",
}
# Collapsed views: promoter_lost stays alone (the one robustly-down class);
# promoter_disconnected trends UP so it is NOT merged into "lost".
GROUP3 = {
    "promoter_lost":              "Promoter\nlost",
    "promoter_disconnected":      "Other\naffected",
    "promoter_proximity_changed": "Other\naffected",
    "context_only":               "Other\naffected",
    "readthrough_exposed":        "Other\naffected",
    "new_promoter_fusion":        "Other\naffected",
    "unaffected":                 "Unaffected",
}
GROUP3_ORDER = ["Promoter\nlost", "Other\naffected", "Unaffected"]
GROUP2 = {c: ("Promoter\nlost" if c == "promoter_lost" else "Promoter\nretained")
          for c in CLASS_ORDER}
GROUP2_ORDER = ["Promoter\nlost", "Promoter\nretained"]

HILITE = "#D55E00"   # Okabe-Ito vermillion: the one robustly-down class
NEUTRAL = "#cfe2f3"  # light steel blue for all other classes


def _stars(p):
    if not np.isfinite(p):
        return ""
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


# ================================================ impact-class plot (former panel d; SI/optional)
def panel_impact(grouping="g3", kind="violin", figsize=(7 / 4, 7 / 4), min_n_star=5,
                 out_name="R5_impact_class_TPM_FC.pdf"):
    """log10 Syn3A/Syn1 transcript fold change by gene_impact_class.

    Final R5 panel c = 3-group violin (Promoter lost / Other affected /
    Unaffected); the full 7-class detail lives in the SI Gene_table.

    grouping : 'g3' (default), 'g2' (Promoter lost / Promoter retained), or
               'full7' (all classes).
    kind     : 'violin' (default) or 'box'.
    Significance stars suppressed for groups with n < min_n_star. Coding only
    (mRNA + pseudo); promoter_lost vs Unaffected stats mirror 09.
    out_name : if None, an exploration name is derived from grouping/kind.
    """
    gt = pd.read_excel(SI_XLSX, sheet_name="Gene_table")
    m = gt[gt["rna_type"].isin(["mRNA", "pseudo"]) &
           gt["gene_impact_class"].notna() &
           gt["TPM_fold_change"].notna() & (gt["TPM_fold_change"] > 0)].copy()
    m["log10FC"] = np.log10(m["TPM_fold_change"])

    if grouping == "g3":
        m["grp"] = m["gene_impact_class"].map(GROUP3)
        order = GROUP3_ORDER; base_label = "Unaffected"
    elif grouping == "g2":
        m["grp"] = m["gene_impact_class"].map(GROUP2)
        order = GROUP2_ORDER; base_label = "Promoter\nretained"
    else:
        m["grp"] = m["gene_impact_class"]
        order = [c for c in CLASS_ORDER if c in set(m["gene_impact_class"])]
        base_label = "unaffected"
        labelmap = CLASS_LABEL

    groups = [g for g in order if g in set(m["grp"])]
    data = [m.loc[m["grp"] == g, "log10FC"].values for g in groups]
    base = m.loc[m["grp"] == base_label, "log10FC"].values

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    pos = range(len(groups))

    if kind == "violin":
        vp = ax.violinplot(data, positions=pos, widths=0.8, showextrema=False)
        for g, body in zip(groups, vp["bodies"]):
            body.set_facecolor(HILITE if "Promoter\nlost" == g else NEUTRAL)
            body.set_alpha(0.85); body.set_edgecolor("0.3"); body.set_linewidth(0.4)
        # thin IQR box + median tick inside each violin
        for i, d in enumerate(data):
            q1, med, q3 = np.percentile(d, [25, 50, 75])
            ax.add_patch(plt.Rectangle((i - 0.04, q1), 0.08, q3 - q1,
                                       facecolor="0.25", edgecolor="none", zorder=4))
            ax.plot([i - 0.13, i + 0.13], [med, med], color="black", lw=1.0, zorder=5)
    else:
        bp = ax.boxplot(data, positions=pos, widths=0.62, showfliers=False,
                        patch_artist=True, medianprops=dict(color="black", linewidth=0.9),
                        boxprops=dict(linewidth=0.6), whiskerprops=dict(linewidth=0.6),
                        capprops=dict(linewidth=0.6))
        for g, patch in zip(groups, bp["boxes"]):
            patch.set_facecolor(HILITE if "Promoter\nlost" == g else NEUTRAL)
            patch.set_alpha(0.9)
        rng = np.random.default_rng(0)
        for i, d in enumerate(data):
            ax.scatter(i + rng.uniform(-0.2, 0.2, size=len(d)), d,
                       s=2.5, color="#525252", alpha=0.35, edgecolors="none", zorder=3)

    ax.axhline(0, color="0.4", linestyle="--", linewidth=0.7, zorder=1)

    log(f"\n[panel d | grouping={grouping} kind={kind}] log10 TPM FC by group")
    ymax = max(np.max(d) for d in data if len(d))
    for i, (g, d) in enumerate(zip(groups, data)):
        if g == base_label or len(d) < 2 or len(base) < 2:
            p = np.nan
        else:
            p = stats.mannwhitneyu(d, base, alternative="two-sided").pvalue
        star = _stars(p) if len(d) >= min_n_star else ""
        med = float(np.median(d)) if len(d) else float("nan")
        log(f"  {g.replace(chr(10),' '):<26} n={len(d):3d}  median FC={10 ** med:.3f}  "
            f"p(vs {base_label.replace(chr(10),' ')})={p:.2e}  star={star or '-'}")
        if star and star != "ns":
            ax.text(i, ymax + 0.12, star, ha="center", va="bottom",
                    fontsize=5, color=(HILITE if "Promoter\nlost" == g else "black"))
        ax.text(i, ax.get_ylim()[0], f"n={len(d)}", ha="center", va="bottom",
                fontsize=4.5, color="0.35")

    xlabels = [labelmap[g] for g in groups] if grouping == "full7" else groups
    ax.set_xticks(list(pos))
    ax.set_xticklabels(xlabels, fontsize=5)
    ax.set_ylabel(r"log$_{10}$ TPM fold change" "\n" r"(Syn3A / Syn1)", fontsize=6)
    ax.tick_params(axis="both", length=2, pad=1.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    if out_name is None:
        tag = {"full7": "full7", "g3": "3grp", "g2": "2grp"}[grouping]
        out_name = f"R5c_TPM_FC_by_impact_class_{tag}_{kind}.pdf"
    out = os.path.join(OUTDIR, out_name)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    log(f"Saved: {out}  (figsize={figsize[0]:.2f} x {figsize[1]:.2f} in)")
    return out


def _hupA_minus10():
    """Quantify the -10 promoter box of the hupA operon (OP_00187, + strand) with
    the same scanner used for the canonical operons (promoter_motif), and log it.
    The box lies in the region deleted in Syn3A (DEL_050), i.e. the lost promoter."""
    syn1op = os.path.join(GR, "..", "Syn1_Operon")
    if syn1op not in sys.path:
        sys.path.insert(0, syn1op)
    import promoter_motif as pm
    ops = pd.read_csv(os.path.join(syn1op, "operons.candidate_blocks.tsv"), sep="\t")
    r = ops[ops["operon_id"] == "OP_00187"].iloc[0]
    tss = int(r["tss"])
    res = pm.scan_minus10(tss, str(r["chrom"]), str(r["strand"]))
    m10_s, m10_e = tss - 12, tss - 7   # -10 hexamer window (+ strand)
    log(f"\n[hupA -10 promoter] operon OP_00187  TSS={tss} ({r['strand']} strand)")
    log(f"  -10 hexamer (TANAAT)   : {res['minus10_6mer']}  "
        f"match={res['match6']} mm={res['mm6']} shift={res['shift6']}")
    log(f"  -10 extended (TNNTANAAT): {res['minus10_9mer']}  "
        f"match={res['match9']} mm={res['mm9']} shift={res['shift9']}")
    log(f"  motif_tier             : {res['motif_tier']}")
    log(f"  -10 window {m10_s}-{m10_e} falls inside DEL_050 (440092-441059) "
        f"-> the promoter was deleted in Syn3A (decapitation)")
    return res


# ---- Syn1 Illumina depth (two-step replicate average, matches avg_sense_TPM) ----
SYN1_ILL_DIR = os.path.join(GR, "..", "Syn1_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph")
SYN1_ILL = [("SRR35996296", 0.25), ("SRR35996297", 0.25), ("SRR35996298", 0.5)]  # mean(96,97) then mean w/ 98
SYN1_CHROM, SYN1_LEN_BP = "CP002027.1", 1_078_809


def _bg_window(path, chrom, win_s, win_e):
    import subprocess
    cov = np.zeros(win_e - win_s)
    out = subprocess.run(["awk", "-F", "\t", f'$1=="{chrom}" && $3>{win_s} && $2<{win_e}', path],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        _, s, e, v = line.split("\t"); s, e, v = int(s), int(e), float(v)
        cov[max(s, win_s) - win_s: min(e, win_e) - win_s] = v
    return cov


def load_syn1_illumina_plus(win_s, win_e):
    cov = np.zeros(win_e - win_s)
    for samp, w in SYN1_ILL:
        cov += w * _bg_window(os.path.join(SYN1_ILL_DIR, f"{samp}.plus.bedGraph"), SYN1_CHROM, win_s, win_e)
    return cov


_SYN1_ILL_MEAN = [None]
def syn1_illumina_mean_total():
    """Genome-wide mean per-base (plus+minus) Syn1 Illumina depth, same replicate weights."""
    if _SYN1_ILL_MEAN[0] is None:
        import subprocess
        tot = 0.0
        for strand in ("plus", "minus"):
            for samp, w in SYN1_ILL:
                p = os.path.join(SYN1_ILL_DIR, f"{samp}.{strand}.bedGraph")
                out = subprocess.run(["awk", "-F", "\t",
                    f'$1=="{SYN1_CHROM}"{{a+=($3-$2)*$4}} END{{print a+0}}', p],
                    capture_output=True, text=True).stdout.strip()
                tot += w * (float(out) if out else 0.0)
        _SYN1_ILL_MEAN[0] = tot / SYN1_LEN_BP
    return _SYN1_ILL_MEAN[0]


def _draw_depth(ax, xg, cov, fill_c, line_c, label, nice_top, logy=False):
    """Depth fill+line on a linear (0..nice_top) or log y-axis. NaN in `cov` marks a
    gap (e.g. a deleted region on the Syn3A track); on a log axis it stays blank, on a
    linear axis it drops to the 0 baseline. A log axis lets a low-expressed operon show
    alongside a high one in the same track."""
    c = np.asarray(cov, dtype=float)
    if logy:
        floor = 0.05
        cc = np.where(np.isfinite(c), np.clip(c, floor, None), np.nan)
        ax.set_yscale("log")
        ax.fill_between(xg, floor, cc, color=fill_c, lw=0, zorder=1)
        ax.plot(xg, cc, color=line_c, lw=0.4, zorder=2)
        mx = float(np.nanmax(cc)) if np.isfinite(cc).any() else 1.0
        top = 10.0 ** np.ceil(np.log10(mx))
        ax.set_ylim(floor, top)
        ticks = [t for t in (0.1, 1, 10, 100) if floor <= t <= top]
        ax.set_yticks(ticks); ax.set_yticklabels([f"{t:g}×" for t in ticks])
    else:
        cc = np.nan_to_num(c, nan=0.0)
        ax.fill_between(xg, 0, cc, color=fill_c, lw=0, zorder=1)
        ax.plot(xg, cc, color=line_c, lw=0.4, zorder=2)
        m = float(cc.max())
        T = nice_top(m) if m > 0 else 1.0
        ax.set_yticks([0, T]); ax.set_yticklabels(["0", f"{T:.0f}×" if T >= 1 else f"{T:g}×"])
        ax.set_ylim(0, T * 1.02)
    ax.set_xlim(0, len(xg))
    ax.set_ylabel(label, fontsize=5, color=line_c)
    ax.tick_params(labelsize=5)
    ax.spines[["top", "right"]].set_visible(False)


# ================================= shared 4-track deletion-junction plotter (panels d, e)
def _junction_panel(win_s, win_e, D1, iso_sel, rel_ticks, fig_w, fig_h,
                    label_override=None, max_iso=8, logy=False):
    """Four tracks on a deletion-junction RELATIVE axis (+ strand; rel = syn1_pos - D1), shared
    by panel d (pdh/acetate operon) and panel e (hupA):
      genes (syn3A deletions shaded across all tracks) | Syn1 PacBio isoforms |
      Syn1 Illumina depth (blue) | Syn3A Illumina depth (red, mapped per-base through the
      retained blocks so deletions read as gaps).
    Isoforms are PacBio (transcript structure); depth is Illumina for BOTH organisms (the
    quantitative standard; PacBio under-samples short genes), each x its own genome mean.
    Returns (fig, axes, cov1, cov3); the caller saves so it can add annotations first."""
    r4dir = os.path.join(GR, "..", "Syn1_Novel_ORF")
    if r4dir not in sys.path:
        sys.path.insert(0, r4dir)
    import R4_track_panels as R4

    win_len = win_e - win_s
    _gl = R4.gene_label
    if label_override is not None:
        R4.gene_label = label_override
    try:
        fig, axes = plt.subplots(4, 1, figsize=(fig_w, fig_h),
                                 height_ratios=[1.0, 0.75, 1.0, 1.0], constrained_layout=True)
        R4.draw_gene_track(axes[0], win_s, win_e, '+')
        R4.draw_isoform_track(axes[1], iso_sel, win_s, win_e, '+', color="#1b6ca8", max_iso=max_iso)
    finally:
        R4.gene_label = _gl

    xg = np.arange(win_len)
    cov1 = load_syn1_illumina_plus(win_s, win_e) / syn1_illumina_mean_total()
    _draw_depth(axes[2], xg, cov1, "#9ecae1", "#3182bd", "Syn1\n(× mean)", R4._nice_top, logy)
    axes[2].set_xticks([])

    # Syn3A depth mapped per-base through the retained blocks; deletions stay NaN (true gaps)
    retained = pd.read_excel(os.path.join(GR, "aln/analysis/genome_reduction_summary.xlsx"))
    retained = retained[retained["Change Case"] == "retained_ordered"]
    cov3 = np.full(win_len, np.nan)
    for _, b in retained.iterrows():
        if pd.isna(b["S2"]):
            continue
        s1, e1, s2 = int(b["S1"]), int(b["E1"]), int(b["S2"])
        lo, hi = max(s1, win_s), min(e1, win_e)
        if lo >= hi:
            continue
        q0 = s2 + (lo - s1)
        cov3[lo - win_s:hi - win_s] = R4.load_syn3a_depth_plus(q0, q0 + (hi - lo))
    cov3 /= R4.syn3a_mean_depth_total()
    _draw_depth(axes[3], xg, cov3, "#f3b0ad", "#c0392b", "Syn3A\n(× mean)", R4._nice_top, logy)
    off = D1 - win_s                                                    # rel -> tx: tx = rel + off
    axes[3].set_xticks([r + off for r in rel_ticks])
    axes[3].set_xticklabels([str(r) for r in rel_ticks], fontsize=5)
    axes[3].set_xlabel('Relative genome position (bp)', fontsize=6)

    for d0, d1 in R4.DELETIONS:                          # deletion shading across all four tracks
        if d1 <= win_s or d0 >= win_e:
            continue
        xa, xb = max(d0, win_s) - win_s, min(d1, win_e) - win_s
        for a in (axes[1], axes[2], axes[3]):
            a.axvspan(xa, xb, facecolor='#e8736a', alpha=0.17, lw=0, zorder=0)
    return fig, axes, cov1, cov3


# ================================================ panel c: two decapitated central-carbon operons
def panel_c(out_name="R5c_central_carbon.pdf", fig_w=7, fig_h=7 / 3):
    """Two adjacent, separately decapitated central-carbon operons on one deletion-junction axis
    (both + strand; rel = syn1_pos - 292905):
      OP_00121 (pdhC/0227-lpdA/0228-pta/0229-ackA/0230): DEL 288391-292905 (4,514 bp) removed the
        promoter (TSS 291897) + 0223-0224 + the PDH E1 subunits pdhA/0225 & pdhB/0226.
      OP_00122 (ptsP/0233-0234-0235; syn3A ptsI/crr): DEL 298422-300803 (2,381 bp) removed the
        promoter (TSS 300106) + 0231 & coaD/0232.
    Both retained gene sets are gene_impact_class promoter_lost and drop in Syn3A (OP_00121 FC
    0.45-0.19; OP_00122 FC 0.41-0.80). PacBio isoforms are filtered to those that SPAN each operon,
    so the two co-transcribed units read as two separate isoform stacks (no isoform crosses between
    them, confirming they are distinct operons). 4-track format; Illumina depth both organisms."""
    r4dir = os.path.join(GR, "..", "Syn1_Novel_ORF")
    if r4dir not in sys.path:
        sys.path.insert(0, r4dir)
    import R4_track_panels as R4
    D1 = 292905
    win_s, win_e = 292300, 303600
    plus = R4.ISO[(R4.ISO.chrom == "CP002027.1") & (R4.ISO.strand == '+')]
    # isoforms that SPAN each operon (co-transcription evidence); top-N of each
    span1 = plus[(plus.start0 <= 293500) & (plus.end0 >= 297500)].sort_values('n_reads', ascending=False).head(8)
    span2 = plus[(plus.start0 <= 300200) & (plus.end0 >= 303200)].sort_values('n_reads', ascending=False).head(8)
    sel = pd.concat([span1, span2], ignore_index=True)
    rel_ticks = [0, 2000, 4000, 6000, 8000, 10000]
    fig, axes, cov1, cov3 = _junction_panel(win_s, win_e, D1, sel, rel_ticks, fig_w, fig_h,
                                            max_iso=len(sel), logy=True)
    out = os.path.join(OUTDIR, out_name)
    fig.savefig(out, dpi=300); plt.close(fig)

    # average x-genome-mean depth per region (retained genes) and per gene, for Illustrator labels
    def _avg(cov, g0, g1):
        return float(np.nanmean(cov[g0 - win_s:g1 - win_s]))
    reg = [("OP_00121 (pdhC-lpdA-pta-ackA)", 292934, 298352),   # retained 0227-0230
           ("OP_00122 (ptsI-crr-0235)",      300901, 303460)]   # retained 0233-0235
    genes_bc = [("pdhC/0227", 292934, 294260), ("lpdA/0228", 294278, 296168),
                ("pta/0229", 296189, 297158), ("ackA/0230", 297170, 298352),
                ("ptsI/0233", 300901, 302623), ("crr/0234", 302704, 303169),
                ("0235", 303169, 303460)]
    log(f"\n[panel c] central-carbon operons OP_00121 + OP_00122 decapitated (DEL 288391-292905 "
        f"& 298422-300803); operon-spanning isoforms {len(span1)}+{len(span2)}; log-y depth -> {out}")
    log("  average depth (x genome mean), retained genes; FC = Syn3A/Syn1:")
    for nm, g0, g1 in reg:
        s1v, s3v = _avg(cov1, g0, g1), _avg(cov3, g0, g1)
        log(f"    REGION {nm:<30s} Syn1 {s1v:6.2f}x  Syn3A {s3v:6.2f}x  FC {s3v/s1v:.3f}")
    for nm, g0, g1 in genes_bc:
        s1v, s3v = _avg(cov1, g0, g1), _avg(cov3, g0, g1)
        log(f"    gene   {nm:<30s} Syn1 {s1v:6.2f}x  Syn3A {s3v:6.2f}x  FC {s3v/s1v:.3f}")
    return out


# ============================================================ panel d: HupA decapitation
def panel_d(out_name="R5d_hupA_operon.pdf", fig_w=7, fig_h=7 / 4):
    """HupA operon (OP_00187, + strand), decapitated: DEL_050 (440092-441059) removed the operon
    promoter (TSS 441031, -10 box 441019-441024) inside gpsA/MMSYN1_0349, so hupA collapses in
    Syn3A. 4-track deletion-junction panel (rel = syn1_pos - deletion_end 441059); isoforms PacBio,
    depth Illumina both organisms; window includes recU/0351 and 0353 for completeness."""
    r4dir = os.path.join(GR, "..", "Syn1_Novel_ORF")
    if r4dir not in sys.path:
        sys.path.insert(0, r4dir)
    import R4_track_panels as R4
    D1 = 441059
    win_s, win_e = 440000, 443000
    hup_s, hup_e = 441113, 441386
    sel = R4.ISO[(R4.ISO.strand == '+') & (R4.ISO.start0 < hup_e) & (R4.ISO.end0 > hup_s) &
                 (R4.ISO.n_reads >= 10)]
    _orig = R4.gene_label
    override = lambda r: 'hupA' if str(r.locus_tag) == 'MMSYN1_0350' else _orig(r)
    rel_ticks = [-1000, -500, 0, 500, 1000, 1500]
    fig, axes, cov1, cov3 = _junction_panel(win_s, win_e, D1, sel, rel_ticks, fig_w, fig_h,
                                            label_override=override)
    hs, he = hup_s - win_s, hup_e - win_s
    out = os.path.join(OUTDIR, out_name)
    fig.savefig(out, dpi=300); plt.close(fig)
    log(f"\n[panel d] HupA operon OP_00187 decapitated (DEL_050; TSS 441031 deleted); depth=Illumina "
        f"both organisms; hupA syn1={np.nanmean(cov1[hs:he]):.1f}x vs syn3A={np.nanmean(cov3[hs:he]):.2f}x mean -> {out}")
    _hupA_minus10()
    return out


# ---- transcript-relative (nt from rpsT/0082 5' end) drawers for the combined b+c panel ----
SYN3A_ISO_TSV = os.path.join(GR, "..", "Syn3A_Transcriptomics/Isoform_Cluster/isoform_clusters_annotated.tsv")
SYN3A_CHROM_BC = "CP016816.2"


def _depth_on_tp(sources, chrom, anchor5p, tp_lo, tp_hi, mean):
    """Per-base depth on the transcript axis (minus strand: genomic = anchor5p - tp).
    `sources` = list of (bedGraph_path, weight) so replicate libraries can be averaged."""
    xg = np.arange(tp_lo, tp_hi)
    g = anchor5p - xg
    g_lo, g_hi = int(g.min()), int(g.max()) + 1
    cov_g = np.zeros(g_hi - g_lo)
    for path, w in sources:
        cov_g += w * _bg_window(path, chrom, g_lo, g_hi)
    return xg, cov_g[g - g_lo] / mean


def _genes_tp(ax, genes, anchor, tp_lo, tp_hi):
    """Minus-strand gene arrows on the transcript axis (5'->3' left->right)."""
    from matplotlib.patches import Polygon
    ax.set_xlim(tp_lo, tp_hi); ax.set_ylim(0, 1.5)
    ax.hlines(0.55, tp_lo, tp_hi, color="black", lw=0.8, zorder=1)
    for g0, g1, name, col in genes:
        xl, xr = anchor - g1, anchor - g0                # minus strand -> arrow points right
        head = min(max(25, (xr - xl) * 0.22), xr - xl)
        v = [(xl, 0.55 - 0.16), (xr - head, 0.55 - 0.16), (xr - head, 0.55 - 0.24), (xr, 0.55),
             (xr - head, 0.55 + 0.24), (xr - head, 0.55 + 0.16), (xl, 0.55 + 0.16)]
        ax.add_patch(Polygon(v, closed=True, facecolor=col, edgecolor="black", lw=0.3, zorder=2))
        ax.text((xl + xr) / 2, 0.55 + 0.30, name, ha="center", va="bottom",
                fontsize=5, color="#333", clip_on=True)
    ax.set_yticks([]); ax.set_xticks([])
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)


def _iso_tp(ax, iso, anchor, tp_lo, tp_hi, nice_top, pack_rows, color, max_iso=8, hl=None):
    """Minus-strand isoform arrows on the transcript axis. The `max_iso` most-abundant
    isoforms form the muted population; `hl` (a DataFrame, e.g. the junction-spanning
    isoform, which may be too rare for the top-N) is force-drawn on top in vermillion
    with a read-count label."""
    from matplotlib.patches import FancyArrowPatch
    iso = iso.sort_values("n_reads", ascending=False).head(max_iso)
    ints, meta = [], []
    for _, r in iso.iterrows():
        a, b = anchor - int(r.end0), anchor - int(r.start0)   # 5'(end0) left, 3'(start0) right
        ints.append((min(a, b), max(a, b))); meta.append(int(r.n_reads))
    rows = pack_rows(ints) if ints else []
    nmax = max(1, max(meta) if meta else 1)
    for (xl, xr), ri, nr in zip(ints, rows, meta):
        lw = float(np.clip(0.3 + 0.7 * np.log10(max(1, nr)), 0.5, 2.6))
        ax.add_patch(FancyArrowPatch((xl, ri), (xr, ri), arrowstyle="-|>", lw=lw, color=color,
                     alpha=min(1.0, 0.5 + 0.5 * nr / nmax), shrinkA=0, shrinkB=0,
                     mutation_scale=4, zorder=2))
    top = max(rows) if rows else 0
    n_hl = 0
    if hl is not None and len(hl):
        for k, (_, r) in enumerate(hl.sort_values("n_reads", ascending=False).iterrows()):
            xl, xr = anchor - int(r.end0), anchor - int(r.start0)
            y = top + 1.3 + k
            ax.add_patch(FancyArrowPatch((xl, y), (xr, y), arrowstyle="-|>", lw=1.7,
                         color="#D55E00", shrinkA=0, shrinkB=0, mutation_scale=5, zorder=5))
            ax.text((xl + xr) / 2, y + 0.5, f"{int(r.n_reads)} reads span 0094 & rpsT/0082",
                    ha="center", va="bottom", fontsize=4.5, color="#D55E00")
            n_hl += 1
        top = top + 1.3 + len(hl)
    ax.set_xlim(tp_lo, tp_hi); ax.set_ylim(-1, top + 1.6)
    ax.set_yticks([]); ax.set_xticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    return n_hl


# ============================================================ panel b: rpsT/0082 partner switch (fusion)
def panel_b(out_name="R5b_rpsT_fusion.pdf", fig_w=14 / 3, fig_h=7 / 3):
    """Combined old panels b + c. rpsT/0082's 5' partner switches from 0083 (Syn1,
    co-transcribed, OP_00043) to 0094 (Syn3A, fused, OP_00050) after DEL_014 (111508-126973,
    15,465 bp; 8 genes 0083-0093) removes the intervening span. Both organisms are aligned on
    the retained rpsT/0082 5' end: x = transcript nt from that 5' end (0082 body positive/right,
    5' partner negative/left; 5'->3' left->right), so 0083 (Syn1) and 0094 (Syn3A) fall in the
    same slot -- the partner swap is a direct top/bottom comparison and the 15.5 kb deletion
    needs no coordinate. Six tracks:
      Syn1 genes (0083 + rpsT/0082) | Syn1 PacBio isoforms | Syn1 Illumina depth (x-mean, blue)
      Syn3A genes (0094 + rpsT/0082) | Syn3A ONT isoforms (junction-spanning in orange) | Syn3A Illumina depth (x-mean, red)
    Depth is Illumina for BOTH organisms (quantitative, shows the rpsT crash relTPM 1.5->0.11);
    PacBio (short-gene under-sampling) and ONT (3'-bias on the 3' gene rpsT) each reverse the
    trend, so they are used only for the isoform/bridging structure, not for depth.
    """
    r4dir = os.path.join(GR, "..", "Syn1_Novel_ORF")
    if r4dir not in sys.path:
        sys.path.insert(0, r4dir)
    import R4_track_panels as R4

    A1, A3 = 111369, 63530             # rpsT/0082 5' end (minus strand, high coord) in Syn1 / Syn3A
    TP_LO, TP_HI = -900, 320
    genes1 = [(111502, 112138, "0083", "#7a7a7a"), (111123, 111369, "rpsT/0082", "#7a7a7a")]
    genes3 = [(63663, 64380, "0094", "#7a7a7a"),   (63284, 63530, "rpsT/0082", "#7a7a7a")]

    i1 = R4.ISO[(R4.ISO.chrom == "CP002027.1") & (R4.ISO.strand == "-") &
                (R4.ISO.start0 < 112138) & (R4.ISO.end0 > 111123) & (R4.ISO.n_reads >= 5)]
    s3iso = pd.read_csv(SYN3A_ISO_TSV, sep="\t")
    i3 = s3iso[(s3iso.chrom == SYN3A_CHROM_BC) & (s3iso.strand == "-") &
               (s3iso.start0 < 64380) & (s3iso.end0 > 63284) & (s3iso.n_reads >= 2)]
    # the lone ONT isoform that spans BOTH gene bodies (0094 64380 & rpsT/0082 63284) = the fusion read
    span3 = i3[(i3.start0 <= 63284) & (i3.end0 >= 64380)]

    fig, axes = plt.subplots(6, 1, figsize=(fig_w, fig_h),
                             height_ratios=[0.55, 1.0, 0.85, 0.55, 1.0, 0.85],
                             constrained_layout=True)
    ag1, ai1, ad1, ag3, ai3, ad3 = axes

    def _yax(ax, cov, color, label):
        m = float(cov.max())
        if m > 0:
            T = R4._nice_top(m)
            ax.set_yticks([0, T]); ax.set_yticklabels(["0", f"{T:.0f}×" if T >= 1 else f"{T:g}×"])
            ax.set_ylim(0, T * 1.02)
        else:
            ax.set_ylim(0, 1)
        ax.set_xlim(TP_LO, TP_HI); ax.set_ylabel(label, fontsize=5, color=color)
        ax.tick_params(labelsize=5); ax.spines[["top", "right"]].set_visible(False)

    # depth = Illumina minus for both organisms (comparable; unbiased for the short 3' rpsT)
    syn1_ill_minus = [(os.path.join(SYN1_ILL_DIR, f"{s}.minus.bedGraph"), w) for s, w in SYN1_ILL]

    # --- Syn1 (PacBio isoforms, Illumina depth) ---
    _genes_tp(ag1, genes1, A1, TP_LO, TP_HI)
    ag1.axvspan(TP_LO, A1 - 111508, facecolor="#e8736a", alpha=0.17, lw=0, zorder=0)  # DEL_014: 0083+ deleted in syn3A
    ag1.text((TP_LO + (A1 - 111508)) / 2, 1.18, "deleted in Syn3A", ha="center", va="top",
             fontsize=4.5, color="#c0392b")
    _iso_tp(ai1, i1, A1, TP_LO, TP_HI, R4._nice_top, R4._pack_rows, "#6baed6", max_iso=8)
    xg1, c1 = _depth_on_tp(syn1_ill_minus, "CP002027.1", A1, TP_LO, TP_HI, syn1_illumina_mean_total())
    ad1.fill_between(xg1, 0, c1, color="#9ecae1", lw=0, zorder=1)
    ad1.plot(xg1, c1, color="#3182bd", lw=0.4, zorder=2)
    _yax(ad1, c1, "#3182bd", "Syn1\n(× mean)"); ad1.set_xticks([])

    # --- Syn3A (ONT isoforms, Illumina depth) ---
    _genes_tp(ag3, genes3, A3, TP_LO, TP_HI)
    n_hl = _iso_tp(ai3, i3, A3, TP_LO, TP_HI, R4._nice_top, R4._pack_rows, "#bdbdbd",
                   max_iso=8, hl=span3)
    xg3, c3 = _depth_on_tp([(R4.SYN3A_DEPTH_MINUS, 1.0)], SYN3A_CHROM_BC, A3, TP_LO, TP_HI,
                           R4.syn3a_mean_depth_total())
    ad3.fill_between(xg3, 0, c3, color="#f3b0ad", lw=0, zorder=1)
    ad3.plot(xg3, c3, color="#c0392b", lw=0.4, zorder=2)
    _yax(ad3, c3, "#c0392b", "Syn3A\n(× mean)")
    ad3.set_xticks([-750, -500, -250, 0, 250])
    ad3.set_xticklabels(["-750", "-500", "-250", "0", "250"], fontsize=5)
    ad3.set_xlabel("Transcript position from rpsT/0082 5′ (nt)", fontsize=6)

    # shared cues: shade the retained rpsT/0082 body across all tracks + 5'-end line;
    # tag each organism block; note the fused partner.
    for ax in axes:
        ax.axvspan(0, 246, facecolor="#dddddd", alpha=0.35, lw=0, zorder=0)
        ax.axvline(0, color="#888888", ls=":", lw=0.6, zorder=0)
    ag1.text(TP_LO, 1.55, "Syn1.0", fontsize=6, fontweight="bold", color="#3182bd", va="bottom")
    ag3.text(TP_LO, 1.55, "Syn3A", fontsize=6, fontweight="bold", color="#c0392b", va="bottom")
    ai1.set_ylabel("PacBio\niso", fontsize=5); ai3.set_ylabel("ONT\niso", fontsize=5)

    out = os.path.join(OUTDIR, out_name)
    fig.savefig(out, dpi=300); plt.close(fig)
    log(f"\n[panel b] rpsT/0082 partner switch 0083(Syn1)->0094(Syn3A) via DEL_014 (15,465 bp); "
        f"Syn1 iso={len(i1)}, Syn3A iso={len(i3)} ({n_hl} junction-spanning shown); "
        f"rpsT depth syn1={float(c1[-TP_LO:-TP_LO+246].mean()):.1f}x vs syn3A={float(c3[-TP_LO:-TP_LO+246].mean()):.2f}x mean -> {out}")
    return out


# ============================================================ panel a: genome-reduction map
def panel_a(out_name="R5a_genome_reduction_map.pdf", figsize=(7 / 3, 7 / 3)):
    """Circular Syn1 -> Syn3A reduction map (matplotlib, born-at-size). Outer ring
    = Syn1 (retained gray, 95 deletions red with arc height ~ log length), inner
    ring = Syn3A (retained gray, novel insertion blue); the one relocated gene lap
    is orange on both rings, linked by a chord. Center = headline count + a
    deletion-length histogram. Reads genome_reduction_summary.xlsx + the BED."""
    import matplotlib.patches as mpatches
    SYN1_LEN, SYN3A_LEN = 1_078_809, 543_379
    df = pd.read_excel(os.path.join(GR, "aln/analysis/genome_reduction_summary.xlsx"))
    retained  = df[df["Change Case"] == "retained_ordered"]
    inserted  = df[df["Change Case"] == "inserted"]
    relocated = df[df["Change Case"] == "retained_relocated"]
    bed = pd.read_csv(os.path.join(GR, "aln/raw/syn1_deleted_regions.bed"),
                      sep="\t", comment="#", header=None)
    dS, dE = bed[1].astype(int).values, bed[2].astype(int).values
    delL = dE - dS
    n_del = len(delL)

    def arc(s, e, total):
        return 2 * np.pi * ((s + e) / 2.0) / total, 2 * np.pi * (e - s) / total

    OUT_BASE, OUT_H = 0.80, 0.055
    IN_BASE, IN_H = 0.46, 0.055
    GRAY, RED, ORANGE, BLUE = "#cfcfcf", "#d62728", "#ff7f0e", "#1f77b4"

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection="polar")   # fill the panel
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_ylim(0, 1.2)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines["polar"].set_visible(False); ax.grid(False)

    # Syn1 outer ring -- retained (gray) then deletions (red, height ~ log length)
    for _, r in retained.iterrows():
        c, w = arc(int(r["S1"]), int(r["E1"]), SYN1_LEN)
        ax.bar(c, OUT_H, width=w, bottom=OUT_BASE, color=GRAY, edgecolor="none", zorder=2)
    lmin, lmax = np.log10(delL.min()), np.log10(delL.max())
    for s, e, L in zip(dS, dE, delL):
        c, w = arc(s, e, SYN1_LEN); w = max(w, np.deg2rad(0.5))
        h = OUT_H + (np.log10(L) - lmin) / (lmax - lmin) * 0.22
        ax.bar(c, h, width=w, bottom=OUT_BASE, color=RED, edgecolor="none", zorder=3)

    # Syn3A inner ring -- retained (gray) then novel insertion (blue, met14p-bearing)
    for _, r in retained.iterrows():
        if pd.isna(r["S2"]) or pd.isna(r["E2"]):
            continue
        c, w = arc(int(r["S2"]), int(r["E2"]), SYN3A_LEN)
        ax.bar(c, IN_H, width=w, bottom=IN_BASE, color=GRAY, edgecolor="none", zorder=2)
    for _, r in inserted.iterrows():
        if int(r["LEN2"]) < 1000:
            continue
        c, w = arc(int(r["S2"]), int(r["E2"]), SYN3A_LEN); w = max(w, np.deg2rad(1.2))
        ax.bar(c, IN_H + 0.05, width=w, bottom=IN_BASE, color=BLUE, edgecolor="none", zorder=3)
        ax.text(c, IN_BASE + IN_H + 0.12, "met14p", ha="center", va="bottom",
                fontsize=5, color=BLUE)

    # lap relocation: orange on both rings + a chord linking the two positions
    rel = relocated.iloc[0]
    c1, w1 = arc(int(rel["S1"]), int(rel["E1"]), SYN1_LEN); w1 = max(w1, np.deg2rad(2.2))
    c2, w2 = arc(int(rel["S2"]), int(rel["E2"]), SYN3A_LEN); w2 = max(w2, np.deg2rad(2.2))
    ax.bar(c1, OUT_H + 0.10, width=w1, bottom=OUT_BASE - 0.02, color=ORANGE,
           edgecolor="black", linewidth=0.3, zorder=5)
    ax.bar(c2, IN_H + 0.10, width=w2, bottom=IN_BASE - 0.02, color=ORANGE,
           edgecolor="black", linewidth=0.3, zorder=5)
    ax.text(c1, OUT_BASE + OUT_H + 0.17, "lap\n(0154)", ha="center", va="bottom",
            fontsize=5, color=ORANGE, linespacing=0.9)

    # sparse Syn1 coordinate ticks
    for kb in range(0, 1000, 250):
        th = 2 * np.pi * (kb * 1000) / SYN1_LEN
        ax.text(th, OUT_BASE + OUT_H + 0.08, f"{kb} kb", ha="center", va="center",
                fontsize=5, color="black", fontweight="bold")

    # center: headline count, centred in the ring
    fig.text(0.5, 0.5, f"{n_del} deletions\n{delL.sum() / 1000:.0f} kb (~50%) removed",
             ha="center", va="center", fontsize=6, fontweight="bold")

    # compact legend (header removed; circle enlarged to fill the panel)
    handles = [mpatches.Patch(color=GRAY, label="retained"),
               mpatches.Patch(color=RED, label="deleted"),
               mpatches.Patch(color=ORANGE, label="relocated"),
               mpatches.Patch(color=BLUE, label="inserted")]
    fig.legend(handles=handles, loc="lower left", fontsize=5, frameon=False,
               handlelength=1.0, borderaxespad=0.3, labelspacing=0.3)

    out = os.path.join(OUTDIR, out_name)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    log(f"\n[panel a] genome-reduction map: {n_del} deletions, "
        f"{delL.sum() / 1000:.0f} kb removed -> {out}")
    return out


# panel_b (syn1 rpsT operon) and panel_c (syn3A fusion) were merged into panel_bc
# (rpsT/0082 partner switch on the shared-0082 transcript axis).

# R5 panel map: a genome map | b rpsT/0082 partner switch (old syn1-b + syn3A-c merged) |
# c two decapitated central-carbon operons (pdh/acetate + PTS) | d hupA decapitation.
# The impact-class violin (former panel d) is kept as panel_impact for the SI / optional.
PANELS = {"a": panel_a, "b": panel_b, "c": panel_c, "d": panel_d, "impact": panel_impact}

if __name__ == "__main__":
    want = sys.argv[1:] or list(PANELS)
    for k in want:
        if k in PANELS:
            PANELS[k]()
        else:
            print(f"[skip] panel '{k}' not implemented yet")
    with open(STATS_TXT, "w") as fh:
        fh.write("\n".join(_log) + "\n")
    print(f"\nStats log -> {STATS_TXT}")
