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


# ============================================================ panel c: impact-class plot
def panel_c(grouping="g3", kind="violin", figsize=(7 / 2, 7 / 6), min_n_star=5,
            out_name="R5c_TPM_FC_by_impact_class.pdf"):
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

    log(f"\n[panel c | grouping={grouping} kind={kind}] log10 TPM FC by group")
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


# ============================================================ panel d: HupA decapitation
def panel_d(out_name="R5d_hupA_operon.pdf", fig_w=7, fig_h=7 / 4):
    """HupA operon (OP_00187, + strand): decapitated because DEL_050 (440092-441059)
    removed its promoter region inside the neighbouring gene gpsA/MMSYN1_0349.
    Reuses the publication single-operon plotter (born-at-size, syn3A-deletion band);
    forces the `hupA` label on MMSYN1_0350 and drops the large flanking pseudogene
    MMSYN1_0354 (and tightens PAD so its deletion band stays out of frame)."""
    syn1op = os.path.join(GR, "..", "Syn1_Operon")
    if syn1op not in sys.path:
        sys.path.insert(0, syn1op)
    import Operon_Visualization as OV

    _gl, _genes, _pad = OV.gene_label, OV.GENES, OV.PAD_BP
    def _force_hupA(r):
        return "hupA" if str(r.get("locus_tag", "")) == "MMSYN1_0350" else _gl(r)
    OV.gene_label = _force_hupA
    OV.GENES = OV.GENES[OV.GENES["locus_tag"] != "MMSYN1_0354"].copy()
    OV.PAD_BP = 60
    try:
        ops = pd.read_csv(os.path.join(syn1op, "operons.candidate_blocks.tsv"), sep="\t")
        row = ops[ops["operon_id"] == "OP_00187"].iloc[0]
        out = os.path.join(OUTDIR, out_name)
        OV.plot_one_operon(row, out, PLOT_DEPTH=True, fig_w=fig_w, fig_h=fig_h)
    finally:
        OV.gene_label, OV.GENES, OV.PAD_BP = _gl, _genes, _pad
    log(f"\n[panel d] HupA operon OP_00187 (decapitation, DEL_050) -> {out}")
    _hupA_minus10()
    return out


def _draw_isoforms_highlight(ax, oc, others, highlight, gap_tx=20):
    """Pack `others` isoforms into greedy rows (muted blue, width ~ log reads),
    then draw `highlight` isoform(s) on top in vermillion with a label. Transcript
    coords via oc.tx_of_genome_pos0 (handles the minus-strand flip)."""
    def _txint(r):
        a = oc.tx_of_genome_pos0(int(r["start0"]))
        b = oc.tx_of_genome_pos0(int(r["end0"]))
        return (min(a, b), max(a, b))
    rows_end = []
    for _, r in others.sort_values("start0").iterrows():
        lo, hi = _txint(r)
        row = next((ri for ri in range(len(rows_end)) if lo > rows_end[ri] + gap_tx), None)
        if row is None:
            rows_end.append(hi); row = len(rows_end) - 1
        else:
            rows_end[row] = hi
        lw = float(np.clip(0.35 + 0.55 * np.log10(r["n_reads"] + 1), 0.35, 1.4))
        ax.plot([lo, hi], [row, row], color="#9ecae1", lw=lw,
                solid_capstyle="round", zorder=2)
    base = max(1, len(rows_end))
    for k, (_, r) in enumerate(highlight.sort_values("start0").iterrows()):
        lo, hi = _txint(r)
        y = base + 0.8 + k
        # keep the thickness ~ abundance convention (this transcript has few reads);
        # it is highlighted by colour + label + position, not by an inflated width
        lw = float(np.clip(0.35 + 0.55 * np.log10(r["n_reads"] + 1), 0.35, 1.4))
        ax.plot([lo, hi], [y, y], color="#D55E00", lw=lw,
                solid_capstyle="round", zorder=6)
        ax.text((lo + hi) / 2, y + 0.55, f"full-span transcript (n={int(r['n_reads'])})",
                ha="center", va="bottom", fontsize=5, color="#D55E00")
    ax.set_ylim(-1, base + 0.8 + len(highlight) + 1.4)
    ax.set_yticks([])


# ============================================================ panel b: fusion DEL_014
def panel_b(out_name="R5b_fusion_DEL014.pdf", figsize=(7 / 2, 7 / 3),
            iso_min_reads=3):
    """Fusion junction DEL_014: in Syn1, OP_00043 (rpsT/0082) and OP_00050 (0094)
    sit ~15.5 kb apart across 8 deleted genes; in Syn3A they are adjacent on the
    minus strand and co-transcribed (37 bridging ONT reads). Shows the Syn3A joined
    locus (0082 and 0094 only; the +strand trmE/0081 is excluded): genes + spanning
    ONT isoforms + depth, with a dashed marker at the new 0094|rpsT junction. Reuses
    the comparison module's Syn3A data and primitives, compact / born-at-size."""
    import Operon_Comparison_Syn1_Syn3A as OC
    OC.GENE_LABEL_FONTSIZE = 6
    OC.LABEL_FONTSIZE = 6
    OC.MAX_ISOFORMS_TO_PLOT = 20

    # dedupe gene labels: "0094/JCVISYN3A_0094" -> "0094"; keep "0082/rpsT"
    def _mk(locus_tag, gene_name):
        gn = (gene_name or "").strip(); lt = (locus_tag or "").strip()
        num = lt.split("_")[-1]
        if gn and not gn.startswith("JCVISYN3A") and gn != lt and gn.lower() != "nan":
            return f"{num}/{gn}"
        return num
    OC.make_gene_label = _mk

    # OP_00043 -> 0082, OP_00050 -> 0094 (the cross-junction pair; 0083/0093 deleted)
    s3_cand = {"JCVISYN3A_0082", "JCVISYN3A_0094"}
    s3g = OC.syn3a_genes[OC.syn3a_genes["locus_tag"].astype(str).isin(s3_cand)].copy()
    s3_chrom = str(s3g["chrom"].iloc[0]); s3_strand = "-"
    s3_s0, s3_e0 = int(s3g["start0"].min()), int(s3g["end0"].max())
    pad = int(OC.PAD_BP_FRAC * (s3_e0 - s3_s0)) + OC.PAD_BP
    ps3, pe3 = s3_s0 - pad, s3_e0 + pad
    oc3 = OC.OperonCoord(chrom=s3_chrom, strand=s3_strand,
                         opid="fusion@syn3A", start0=s3_s0, end0=s3_e0)
    # isoforms over the locus: show the population (muted) and HIGHLIGHT the
    # transcript(s) spanning the entire body of BOTH genes (rpsT/0082 63285-63530
    # and 0094 63664-64380) -- the n_span fusion cluster (ISO_075840).
    G_LO, G_HI = 63285, 64380
    over = ((OC.syn3a_isoforms["chrom"].astype(str) == s3_chrom) &
            (OC.syn3a_isoforms["strand"].astype(str) == s3_strand) &
            (OC.syn3a_isoforms["start0"] < s3_e0) & (OC.syn3a_isoforms["end0"] > s3_s0))
    allf = OC.syn3a_isoforms[over].copy()
    full = allf[(allf["start0"] <= G_LO) & (allf["end0"] >= G_HI)].copy()
    others = allf[(~allf["isoform_id"].isin(full["isoform_id"])) &
                  (allf["n_reads"] >= iso_min_reads)].copy()
    others = others.sort_values("n_reads", ascending=False).head(30)
    d3 = OC.subset_intervals(OC.syn3a_depth_minus, s3_chrom, ps3, pe3)

    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 2.4, 1.0])
    axg = fig.add_subplot(gs[0])
    axi = fig.add_subplot(gs[1], sharex=axg)
    axd = fig.add_subplot(gs[2], sharex=axg)

    # draw genes; draw_gene_arrows adds a Syn3A genomic-position secondary axis
    # below the arrows -- keep it (relabel + shrink fonts to OUTPUT.md)
    pre = set(fig.axes)
    OC.draw_gene_arrows(axg, oc3, s3g)
    gaxes = [x for x in fig.axes if x not in pre]          # the genome twiny axis
    for gax in gaxes:
        # put the Syn3A genome axis at the TOP of the gene track (as in panel d)
        gax.spines["bottom"].set_visible(False)
        gax.spines["top"].set_visible(True)
        gax.spines["top"].set_position(("outward", 2))
        gax.xaxis.set_ticks_position("top")
        gax.xaxis.set_label_position("top")
        gax.tick_params(axis="x", which="both", top=True, labeltop=True,
                        bottom=False, labelbottom=False, labelsize=5)
        gax.set_xlabel("Syn3A genome position (bp)", fontsize=6)
    axg.spines["top"].set_visible(False)   # drop the gene panel's own top spine
    _draw_isoforms_highlight(axi, oc3, others, full)
    OC.draw_depth(axd, oc3, d3, ps3, pe3, s3_strand)
    log(f"\n[panel b] isoforms shown: {len(others)} others + {len(full)} full-span "
        f"highlighted (>= {iso_min_reads} reads for others)")

    # headroom so the gene labels + deletion note (above the arrows) are not clipped
    yl = axg.get_ylim()
    axg.set_ylim(yl[0], yl[1] + 0.9 * (yl[1] - yl[0]))

    # deletion scar at the new 0094|rpsT junction: the 8 genes (~15.5 kb) between
    # them in Syn1 are absent from Syn3A (no coordinate here), so mark + label the
    # excised block rather than drawing absent genes.
    jx = oc3.tx_of_genome_pos0((63530 + 63664) // 2)
    for a in (axg, axi, axd):
        a.axvline(jx, color="#c0392b", linestyle="--", linewidth=0.7, zorder=6)
    axg.text(jx, axg.get_ylim()[1] * 0.99, "8 genes, ~15.5 kb\ndeleted",
             ha="center", va="top", fontsize=5, color="#c0392b", linespacing=0.9)

    x_left = oc3.tx_of_genome_pos0(pe3)   # minus strand: high genome -> tx 0
    x_right = oc3.tx_of_genome_pos0(ps3)
    lo, hi = min(x_left, x_right), max(x_left, x_right)
    for a in [axg, axi, axd] + gaxes:     # incl. genome twiny so it stays aligned
        a.set_xlim(lo, hi)
    # transcript ticks only on the bottom (depth) axis; the Syn3A genome axis sits
    # under the gene arrows, so axg/axi carry no transcript ticks
    for a in (axg, axi):
        a.tick_params(axis="x", which="both", bottom=False, top=False,
                      labelbottom=False, labeltop=False)
        a.set_xlabel("")
    axd.tick_params(axis="x", which="both", top=False, labeltop=False)
    axd.set_xlabel("Transcript coordinate (nt)", fontsize=6)
    axd.set_ylabel("ONT depth", fontsize=6)
    axi.set_ylabel("ONT isoforms", fontsize=6)
    for a in (axi, axd):
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)

    out = os.path.join(OUTDIR, out_name)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    log(f"[panel b] fusion DEL_014 (Syn3A joined OP_00043|OP_00050) -> {out}")
    return out


# ============================================================ panel a: genome-reduction map
def panel_a(out_name="R5a_genome_reduction_map.pdf", figsize=(7 / 2, 7 / 2)):
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

    # center: headline count + deletion-length histogram
    fig.text(0.5, 0.62, f"{n_del} deletions\n{delL.sum() / 1000:.0f} kb (~50%) removed",
             ha="center", va="center", fontsize=6, fontweight="bold")
    axh = fig.add_axes([0.405, 0.40, 0.19, 0.115])
    edges = [50, 200, 1000, 5000, 20000, int(delL.max()) + 1]
    labels = ["<0.2", "0.2-1", "1-5", "5-20", ">20"]
    counts = [int(((delL >= edges[i]) & (delL < edges[i + 1])).sum()) for i in range(5)]
    axh.bar(range(5), counts, color=RED, width=0.82)
    axh.set_xticks(range(5)); axh.set_xticklabels(labels, fontsize=4)
    axh.set_yticks([0, max(counts)]); axh.tick_params(labelsize=4, length=1.5, pad=1)
    axh.set_xlabel("deletion size (kb)", fontsize=4.5, labelpad=1)
    for sp in ("top", "right"):
        axh.spines[sp].set_visible(False)
    axh.spines["left"].set_linewidth(0.4); axh.spines["bottom"].set_linewidth(0.4)

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
        f"{delL.sum() / 1000:.0f} kb; length bins {counts} -> {out}")
    return out


PANELS = {"a": panel_a, "b": panel_b, "c": panel_c, "d": panel_d}

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
