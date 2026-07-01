#!/usr/bin/env python
"""
R6_figure_panels.py  --  born-at-size panels for Results figure R6
("Minimization reallocates transcription toward translation").

Panel inventory (R6 caption):
  a  mRNA-pool composition by Secondary function  -> reuse 09's
        Compare_RNA_Protein/mRNA_pool_composition_by_secondary.pdf  (already born-at-size)
  b  mRNA fold-change vs absolute change, r-proteins highlighted green, with a
     Syn1-vs-Syn3A relative-mRNA correlation inset (top-left) showing the hierarchy
     is conserved (r=0.84) yet most genes shifted down  [THIS SCRIPT]
     (replaces the old tertiary-share dumbbell, per the 2026-07 R6 restructure)
  c  mRNA + protein fold change of RNA polymerase, the degradosome, and the
     central-carbon enzymes, transposed to a landscape (7 x 7/6) lollipop  [THIS SCRIPT]
  d  The 21-gene ~11 kb r-protein operon + its swapped upstream neighbour, MERGED across
     both cells on one transcript axis anchored on the shared rpsJ/0672 5' end (Syn1 genes+
     Illumina depth over Syn3A genes+Illumina depth). Replaces old d (operon structure) and
     old e (tRNA junction, R6_panel_e_trna_rptn.py -> superseded).  [THIS SCRIPT]

Panels b,c read syn1_vs_syn3a_RNA_protein.tsv (fold changes) +
macromolecule_complex_abundance.tsv (RNAP, degradosome limiting-subunit estimates);
panel d reads the Syn1/Syn3A Illumina minus-strand bedGraphs (via R5's _depth_on_tp helpers
and R4's Syn3A depth) + the Syn1 GFF; recompute nothing.
"""
import os
import sys
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None
from scipy import stats

mpl.rcParams.update({
    "font.size": 7, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
    "ytick.labelsize": 6, "legend.fontsize": 6, "pdf.fonttype": 42, "ps.fonttype": 42,
    # mathtext in Arial (custom fontset) so $...$ labels embed ArialMT, not DejaVu/Cmsy10
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "mathtext.default": "rm",
})

GR     = os.path.dirname(os.path.abspath(__file__))
CRP    = os.path.join(GR, "Compare_RNA_Protein")
OUTDIR = os.path.join(GR, "R6_panels")
os.makedirs(OUTDIR, exist_ok=True)

TPM_C, IPM_C = "#3182bd", "#e6550d"   # transcript = blue, protein = orange
GREEN = (16 / 255, 130 / 255, 60 / 255)   # ribosomal-protein highlight


# ===================================== panel b: mRNA fold-change vs absolute change,
#                                       r-proteins highlighted, with a syn1-vs-syn3A
#                                       correlation inset (the "most genes shift down" cue)
def panel_b(out_name="R6b_mRNA_FC_vs_absChange.pdf", figsize=(14 / 3, 7 / 2)):
    """mRNA fold change (x, log) vs absolute change (y, symlog) for the retained
    coding pool. Base dots are black with opacity ramped by the Syn1 baseline
    (abundant genes = the ones that dominate absolute change stand out); the
    ribosomal proteins are recoloured green to show translation is the up-mover.
    An inset (top-left) is the Syn1-vs-Syn3A relative-mRNA log-log correlation with
    the y=x diagonal: the hierarchy is conserved (r=0.84) yet most genes sit below
    the line, so the retained pool shifted down to fund the r-protein rise."""
    df = pd.read_csv(os.path.join(CRP, "syn1_vs_syn3a_RNA_protein.tsv"), sep="\t")
    df = df[df["rna_type"] == "mRNA"].copy()
    df["locus_num"] = df["locus_syn1"].str.extract(r"(\d+)$").astype(int)   # trailing digits (avoid the '1' in MMSYN1)
    gn = df["gene_name"].fillna("").astype(str)
    gp = df["gene_product"].fillna("").astype(str).str.lower()
    df["is_rp"] = gn.str.match(r"rp[slm]") | gp.str.contains("ribosomal protein")

    sub = df[df["TPM_fold_change"].notna() & (df["TPM_fold_change"] > 0) &
             df["TPM_abs_change"].notna() &
             df["relTPM_syn1"].notna() & (df["relTPM_syn1"] > 0)].copy()
    x, y = sub["TPM_fold_change"], sub["TPM_abs_change"]
    nz = y[y != 0].abs()
    linthresh = max(1e-3, float(np.nanmedian(nz))) if len(nz) else 1e-3

    # base = black, opacity ramped by log10 Syn1 baseline; r-proteins = solid green.
    logb = np.log10(sub["relTPM_syn1"])
    lo, hi = float(logb.min()), float(logb.max())
    norm = (logb - lo) / (hi - lo) if hi > lo else pd.Series(0.5, index=logb.index)
    rgba = np.zeros((len(sub), 4))
    rgba[:, 3] = 0.12 + 0.78 * norm.to_numpy()          # black base, alpha ramp
    rpmask = sub["is_rp"].fillna(False).to_numpy().astype(bool)
    rgba[rpmask, 0], rgba[rpmask, 1], rgba[rpmask, 2] = GREEN
    rgba[rpmask, 3] = 0.85                               # keep the 52 r-proteins visible

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.scatter(x, y, s=10, facecolors=rgba, edgecolors="none", zorder=2)
    ax.axvline(1, color="black", ls="--", lw=0.7, zorder=1)
    ax.axhline(0, color="black", ls="--", lw=0.7, zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=linthresh)
    # exact axis labels from the old TPM_FC_vs_absChange.pdf (mathtext now in Arial)
    ax.set_xlabel(r"TPM fold change  ($\mathrm{rel}_{\mathrm{syn3A}}/\mathrm{rel}_{\mathrm{syn1}}$)", fontsize=7)
    ax.set_ylabel(r"TPM absolute change  ($\mathrm{rel}_{\mathrm{syn3A}}-\mathrm{rel}_{\mathrm{syn1}}$)", fontsize=7)
    mlog = float(np.log10(x).abs().max()) * 1.05
    ax.set_xlim(10 ** (-mlog), 10 ** mlog)
    ymax = float(y.abs().max()) * 1.15
    ax.set_ylim(-ymax, ymax)
    ax.tick_params(axis="both", length=2, pad=1.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(handles=[Line2D([0], [0], marker="o", ls="", color=GREEN, ms=4,
                              label=f"ribosomal protein (n={int(rpmask.sum())})")],
              loc="lower right", fontsize=5, frameon=False, handletextpad=0.2)

    # colorbar (right): black alpha-ramp mirroring the dot opacity = log10 Syn1 baseline
    cmap = LinearSegmentedColormap.from_list("black_alpha", [(0, 0, 0, 0.12), (0, 0, 0, 0.90)])
    sm = plt.cm.ScalarMappable(norm=Normalize(vmin=lo, vmax=hi), cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label("log10 Syn1 baseline (rel-TPM)", fontsize=6)
    cbar.ax.tick_params(labelsize=5, length=2)

    # annotate exactly as the old pdf: top-5 absolute movers (both signs) + top-5 fold-change
    # outliers (both extremes), union de-duplicated; label = 4-letter gene name, else 4-digit locus.
    top = pd.concat([sub[sub["TPM_abs_change"] > 0].nlargest(5, "TPM_abs_change"),
                     sub[sub["TPM_abs_change"] < 0].nsmallest(5, "TPM_abs_change"),
                     sub.nlargest(5, "TPM_fold_change"),
                     sub.nsmallest(5, "TPM_fold_change")]).drop_duplicates(subset="locus_num")

    def _lab(r):
        g = str(r["gene_name"]).strip()
        if g and g.lower() != "nan":
            return g.replace("/", " ").replace(";", " ").split()[0]   # primary symbol, e.g. folE/yqfO-like -> folE
        return f"{int(r['locus_num']):04d}"
    texts = [ax.text(r["TPM_fold_change"], r["TPM_abs_change"], _lab(r), fontsize=5)
             for _, r in top.iterrows()]
    if adjust_text is not None:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.4))

    # ---- inset (top-left): Syn1 vs Syn3A relative-mRNA log-log correlation ----
    both = df[(df["relTPM_syn1"] > 0) & (df["relTPM_syn3a"] > 0)]
    lx, ly = np.log10(both["relTPM_syn1"]), np.log10(both["relTPM_syn3a"])
    r, _ = stats.pearsonr(lx, ly)
    below = float((both["relTPM_syn3a"] < both["relTPM_syn1"]).mean() * 100)
    axi = ax.inset_axes([0.135, 0.60, 0.34, 0.38])   # shifted right so its ylabel clears the outer y-axis
    axi.scatter(lx, ly, s=1.6, color="0.35", alpha=0.35, edgecolors="none", zorder=2)
    lim = [min(lx.min(), ly.min()) - 0.1, max(lx.max(), ly.max()) + 0.1]
    axi.plot(lim, lim, color="black", ls="--", lw=0.6, zorder=3)   # y = x
    axi.set_xlim(lim); axi.set_ylim(lim)
    axi.set_xlabel(r"rel. TPM$_{\mathrm{syn1}}$", fontsize=5, labelpad=1)
    axi.set_ylabel(r"rel. TPM$_{\mathrm{syn3A}}$", fontsize=5, labelpad=1)
    axi.tick_params(axis="both", length=1.5, labelsize=4, pad=1)
    axi.text(0.05, 0.95, f"r = {r:.2f}\n{below:.0f}% below syn1",
             transform=axi.transAxes, va="top", ha="left", fontsize=4.5)
    for sp in ("top", "right"):
        axi.spines[sp].set_visible(False)

    out = os.path.join(OUTDIR, out_name)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[panel b] mRNA FC-vs-absChange + corr inset (r={r:.3f}, {below:.0f}% below) -> {out}")
    return out


# ===================================== panel c: complex + enzyme fold-change lollipop
def panel_c(out_name="R6c_complex_enzyme_FC.pdf", figsize=(7, 7 / 6)):
    """Vertical lollipop (landscape): per entity, mRNA (TPM) and protein (iPM) fold
    change (Syn3A/Syn1) on a log y-axis with a reference line at 1. RNA polymerase
    and the degradosome (limiting-subunit estimates) on the left, then the
    central-carbon enzymes left-to-right in pathway order."""
    cx = pd.read_csv(os.path.join(CRP, "macromolecule_complex_abundance.tsv"), sep="\t")
    df = pd.read_csv(os.path.join(CRP, "syn1_vs_syn3a_RNA_protein.tsv"), sep="\t")

    def cxrow(name):
        r = cx[cx["complex"] == name].iloc[0]
        return float(r["TPM_fold_change"]), float(r["iPM_fold_change"])

    PRIM = {"Genetic Information Processing": "#3b6db3", "Metabolism": "#3f9e5a"}  # panel-b palette
    RELABEL = {"0607": "GapA", "0451": "GapN"}   # the two G3PDHs: GapA (phosphorylating) vs GapN (non-phos)

    def loc_row(loc):
        r = df[df["locus_syn1"] == f"MMSYN1_{loc}"].iloc[0]
        g = RELABEL.get(loc, str(r["gene_name"]))
        disp = g[0].upper() + g[1:]                            # capitalized = enzyme (protein) name
        return (f"{disp}\n{loc}",                              # name over locus (2nd row)
                float(r["TPM_fold_change"]), float(r["iPM_fold_change"]))

    # ordered groups top -> bottom; metabolic block runs sugar-in to ATP products.
    # header = tertiary-function annotation; label colour = Primary family (as in panel b).
    GIP, MET = "Genetic Information Processing", "Metabolism"
    PTS = ["0233", "0694", "0234", "0779"]                    # ptsI, ptsH, crr, ptsG
    GLY = ["0445", "0220", "0131", "0727", "0607", "0451",    # pgi pfkA fbaA tpiA GapA GapN
           "0606", "0729", "0213", "0221"]                    # pgk pgm eno pyk   (pathway order)
    PYR = ["0227", "0229", "0230", "0475"]                    # pdhC pta ackA ldh
    machines = [("RNAP",) + cxrow("RNA polymerase"),
                ("Degradosome",) + cxrow("Degradosome")]
    groups = [("Gene Expression",        GIP, machines),
              ("Carbohydrate transport", MET, [loc_row(l) for l in PTS]),
              ("Glycolysis",             MET, [loc_row(l) for l in GLY]),
              ("Pyruvate metabolism",    MET, [loc_row(l) for l in PYR])]

    # abbreviated group headers so they fit above a landscape group
    HDR = {"Gene Expression": "Gene expr.", "Carbohydrate transport": "Carb. transport",
           "Glycolysis": "Glycolysis", "Pyruvate metabolism": "Pyruvate metab."}
    GAP = 0.9
    xs, labels, lab_colors, data, headers = [], [], [], [], []
    x = 0.0
    for hdr, fam, ents in groups:
        col = PRIM[fam]
        x0 = x
        for lab, tpm, ipm in ents:
            xs.append(x); labels.append(lab.replace("\n", " ")); lab_colors.append(col)
            data.append((tpm, ipm)); x += 1.0
        headers.append(((x0 + x - 1) / 2.0, HDR.get(hdr, hdr), col))   # group centre
        x += GAP

    fig, ax = plt.subplots(figsize=figsize)
    ax.axhline(1.0, color="0.45", ls="--", lw=0.7, zorder=1)
    for xx, (tpm, ipm) in zip(xs, data):
        ax.plot([xx, xx], [tpm, ipm], color="0.75", lw=0.8, zorder=2)
        ax.scatter([xx], [tpm], facecolor=TPM_C, edgecolor="none", s=11, zorder=3)
        ax.scatter([xx], [ipm], facecolor=IPM_C, edgecolor="none", s=11, zorder=3)
    # tertiary-function header, centred above each group, family-coloured
    for hx, hdr, col in headers:
        ax.text(hx, 1.015, hdr, transform=ax.get_xaxis_transform(), ha="center",
                va="bottom", fontsize=5, fontstyle="italic", fontweight="bold", color=col)

    ax.set_yscale("log")
    ax.set_ylim(0.15, 2.2)
    ax.set_yticks([0.25, 0.5, 1, 2])
    ax.set_yticklabels(["0.25", "0.5", "1", "2"])
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=5, rotation=45, ha="right", va="top", rotation_mode="anchor")
    for tick, col in zip(ax.get_xticklabels(), lab_colors):  # colour labels by Primary family
        tick.set_color(col)
    ax.set_xlim(min(xs) - 0.7, max(xs) + 0.7)
    ax.set_ylabel("Fold Change\n(Syn3A/Syn1)", fontsize=6)
    ax.tick_params(axis="both", length=2, pad=1.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    handles = [Line2D([0], [0], marker="o", color=TPM_C, ls="", ms=4, label="mRNA"),
               Line2D([0], [0], marker="o", color=IPM_C, ls="", ms=4, label="Protein")]
    ax.legend(handles=handles, fontsize=5, loc="upper right", frameon=False, ncol=2,
              handletextpad=0.2, columnspacing=0.8)

    # explicit margins: wider left gutter for the y-axis label/ticks, room below for 45deg labels
    fig.subplots_adjust(left=0.085, right=0.995, top=0.83, bottom=0.34)

    out = os.path.join(OUTDIR, out_name)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[panel c] complex+enzyme FC lollipop -> {out}")
    return out


# ===================================== panel d: the 21-gene rPtn operon + its swapped
#                                       upstream neighbour (merged old d + e)
def _load_operon_genes():
    """rPtn operon 0652-0672 from the Syn1 GFF -> [(start1, end1, gene), ...] (1-based)."""
    gff = os.path.join(GR, "..", "Genomes_Input", "syn1.genes.gff3")
    out = []
    for ln in open(gff):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.search(r"locus_tag=MMSYN1_(\d+)", f[8])
        if not m:
            continue
        n = int(m.group(1))
        if 652 <= n <= 672:
            gm = re.search(r"gene=([^;]+)", f[8])
            out.append((int(f[3]), int(f[4]), gm.group(1) if gm else f"{n:04d}"))
    return out


def panel_d(out_name="R6d_rPtn_operon.pdf", figsize=(7, 7 / 3)):
    """The ~11 kb 21-gene ribosomal-protein operon (rpsJ/0672 -> secY/0652, minus strand)
    with its upstream neighbour, merged across both cells on ONE transcript axis anchored
    on the shared operon 5' end / TSS (rel 0; operon body positive/right, upstream negative/left,
    5'->3' left->right; rpsJ starts at ~+77 after the 5' UTR). Four tracks, top->down:
      Syn1 genes  (upstream = dhaK/0673, co-directional)      | Syn1 Illumina depth (x mean)
      Syn3A genes (upstream = the relocated 4-tRNA operon)     | Syn3A Illumina depth (x mean)
    The retained operon aligns 1:1 between cells; only the upstream neighbour changed
    (dhaK -> tRNAs, pulled in by the flanking deletions), and no read-through bridges the
    silent inter-operon gap in Syn3A. Depth is Illumina (quantitative) for both, x genome-mean.
    Exports the mean normalized operon-body depth (rel 0..10853) to R6de_rPtn_operon_depth.txt."""
    import R5_figure_panels as R5
    r4dir = os.path.join(GR, "..", "Syn1_Novel_ORF")
    if r4dir not in sys.path:
        sys.path.insert(0, r4dir)
    import R4_track_panels as R4
    syn1op = os.path.join(GR, "..", "Syn1_Operon")
    if syn1op not in sys.path:
        sys.path.insert(0, syn1op)
    import promoter_motif as pm

    TSS = 806176                      # OP_00341 TSS (Syn1, minus strand) = the operon 5' end / origin
    A1, A3 = TSS, TSS - 386382        # anchor = operon 5' (TSS); Syn3A TSS via the retained-block offset (419794)
    OFF = A1 - A3                      # 386382; Syn3A operon coord = Syn1 coord - OFF (retained, identical)
    TP_LO, TP_HI = -2000, 11100
    OP_LO, OP_HI = 0, A1 - 795222     # operon body on the transcript axis: TSS (0) -> operon 3' end / tts (10954)

    OPC, GRAY, TEAL = "#8aa9c8", "#b0b0b0", "#2c9e8f"
    og = _load_operon_genes()
    # arrows only (names blank); gene names are drawn separately, alternating above/below the axis
    genes1 = [(s, e, "", OPC) for s, e, nm in og] + [(806351, 807349, "", GRAY)]
    genes3 = [(s - OFF, e - OFF, "", OPC) for s, e, nm in og]
    TRNAS = [(420560, 420635), (420648, 420723), (420731, 420806), (420814, 420889)]  # Thr Val Glu Asn
    genes3 += [(s, e, "", TEAL) for s, e in TRNAS]
    lab1 = [(s, e, nm) for s, e, nm in og] + [(806351, 807349, "dhaK")]   # operon names + dhaK (Syn1)
    lab3 = [(s - OFF, e - OFF, nm) for s, e, nm in og]                    # same operon names (Syn3A)

    def _alt_labels(ax, items, anchor):
        """Gene names centred on each arrow, alternating above/below the axis so the 21
        tightly-packed operon genes stay legible (left->right order sets the alternation)."""
        for i, (g0, g1, nm) in enumerate(sorted(items, key=lambda t: anchor - (t[0] + t[1]) / 2)):
            if not nm:
                continue
            xc = anchor - (g0 + g1) / 2
            if i % 2 == 0:
                ax.text(xc, 0.92, nm, ha="center", va="bottom", fontsize=4.2, color="#333", clip_on=False)
            else:
                ax.text(xc, 0.18, nm, ha="center", va="top", fontsize=4.2, color="#333", clip_on=False)

    # mRNA-pool share of the operon (coding pool), both cells
    dfp = pd.read_csv(os.path.join(CRP, "syn1_vs_syn3a_RNA_protein.tsv"), sep="\t")
    opm = dfp[dfp["locus_syn1"].isin([f"MMSYN1_{n:04d}" for n in range(652, 673)])]
    mrna = dfp[dfp["rna_type"] == "mRNA"]
    sh1 = 100 * opm["relTPM_syn1"].sum() / mrna["relTPM_syn1"].sum()
    sh3 = 100 * opm["relTPM_syn3a"].sum() / mrna["relTPM_syn3a"].sum()
    m10 = pm.scan_minus10(TSS, "CP002027.1", "-")   # operon -10 box (retained in both cells)

    fig, (ag1, ad1, ag3, ad3) = plt.subplots(
        4, 1, figsize=figsize, height_ratios=[0.6, 1.35, 0.6, 1.35], constrained_layout=True)

    syn1_ill_minus = [(os.path.join(R5.SYN1_ILL_DIR, f"{s}.minus.bedGraph"), w) for s, w in R5.SYN1_ILL]

    def _depthrow(ax, xg, cov, fill, line, label):
        ax.fill_between(xg, 0, cov, color=fill, lw=0, zorder=1)
        ax.plot(xg, cov, color=line, lw=0.4, zorder=2)
        m = float(np.nanmax(cov)) if len(cov) else 0.0
        T = R4._nice_top(m) if m > 0 else 1
        ax.set_ylim(0, T * 1.03); ax.set_yticks([0, T])
        ax.set_yticklabels(["0", f"{T:.0f}×" if T >= 1 else f"{T:g}×"], fontsize=5)
        ax.set_xlim(TP_LO, TP_HI); ax.set_ylabel(label, fontsize=5, color=line)
        ax.tick_params(labelsize=5, length=2, pad=1)
        ax.spines[["top", "right"]].set_visible(False)

    # --- Syn1 ---
    R5._genes_tp(ag1, genes1, A1, TP_LO, TP_HI)
    _alt_labels(ag1, lab1, A1)
    ag1.axvspan(TP_LO, A1 - 806355, facecolor="#e8736a", alpha=0.15, lw=0, zorder=0)  # DEL_074 (dhaK+) deleted in Syn3A; starts at -179
    ag1.text(TP_LO + 40, 0.06, "deleted in Syn3A", ha="left", va="bottom", fontsize=4.5, color="#c0392b")
    xg1, c1 = R5._depth_on_tp(syn1_ill_minus, R5.SYN1_CHROM, A1, TP_LO, TP_HI, R5.syn1_illumina_mean_total())
    _depthrow(ad1, xg1, c1, "#9ecae1", "#3182bd", "Syn1\n(× mean)")

    # --- Syn3A ---
    R5._genes_tp(ag3, genes3, A3, TP_LO, TP_HI)
    _alt_labels(ag3, lab3, A3)
    ag3.text(-930, 1.02, "tRNAs\n(Thr Val Glu Asn)", ha="center", va="bottom", fontsize=4.5,
             color=TEAL, linespacing=0.9)
    xg3, c3 = R5._depth_on_tp([(R4.SYN3A_DEPTH_MINUS, 1.0)], "CP016816.2", A3, TP_LO, TP_HI,
                              R4.syn3a_mean_depth_total())
    _depthrow(ad3, xg3, c3, "#f3b0ad", "#c0392b", "Syn3A\n(× mean)")
    ad3.axvspan(A3 - 420350, min(0, A3 - 419784), facecolor="#c0392b", alpha=0.08, lw=0, zorder=0)  # silent inter-operon gap

    # shared cue: 5'-end line across all tracks; organism tags
    for ax in (ag1, ad1, ag3, ad3):
        ax.axvline(0, color="#888888", ls=":", lw=0.6, zorder=0)
    ag1.text(TP_LO, 1.55, "Syn1.0", fontsize=6, fontweight="bold", color="#3182bd", va="bottom")
    ag3.text(TP_LO, 1.55, "Syn3A", fontsize=6, fontweight="bold", color="#c0392b", va="bottom")

    for ax in (ag1, ad1, ag3):
        ax.tick_params(axis="x", labelbottom=False, bottom=False)
    ad3.set_xticks([-1000, 0, 2000, 4000, 6000, 8000, 10000])
    ad3.set_xticklabels(["−1000", "0", "2000", "4000", "6000", "8000", "10000"], fontsize=5)
    ad3.set_xlabel("Relative transcript position (nt)", fontsize=6)

    out = os.path.join(OUTDIR, out_name)
    fig.savefig(out, dpi=300)
    plt.close(fig)

    # export mean normalized operon-body depth (rel 0..10853)
    mask1 = (xg1 >= OP_LO) & (xg1 < OP_HI)
    mask3 = (xg3 >= OP_LO) & (xg3 < OP_HI)
    s1v, s3v = float(np.nanmean(c1[mask1])), float(np.nanmean(c3[mask3]))
    depth_txt = os.path.join(OUTDIR, "R6de_rPtn_operon_depth.txt")
    with open(depth_txt, "w") as fh:
        fh.write("21-gene rPtn operon OP_00341 (rpsJ/0672 -> secY/0652, minus strand)\n")
        fh.write("=" * 64 + "\n")
        fh.write(f"mean Illumina depth over the operon body (rel 0..{int(OP_HI)} nt, TSS -> tts), "
                 "normalized to genome-mean coverage (x mean):\n")
        fh.write(f"  Syn1  : {s1v:.3f}x\n  Syn3A : {s3v:.3f}x\n  fold change (Syn3A/Syn1): {s3v / s1v:.3f}\n\n")
        fh.write("coding mRNA-pool share (sum relTPM of the 21 genes / total coding relTPM):\n")
        fh.write(f"  Syn1  : {sh1:.2f}%\n  Syn3A : {sh3:.2f}%   (share fold change {sh3 / sh1:.2f})\n\n")
        fh.write(f"promoter -10 box (OP_00341 TSS {TSS}, minus strand; retained in both cells, "
                 f"nearest deletion DEL_074 starts 179 bp upstream):\n")
        fh.write(f"  -10 hexamer (TANAAT)    : {m10['minus10_6mer']}  "
                 f"match={m10['match6']} mm={m10['mm6']} shift={m10['shift6']}\n")
        fh.write(f"  -10 extended (TNNTANAAT): {m10['minus10_9mer']}  "
                 f"match={m10['match9']} mm={m10['mm9']} shift={m10['shift9']}\n")
        fh.write(f"  motif_tier              : {m10['motif_tier']}\n")
    print(f"[panel d] merged rPtn-operon panel -> {out}")
    print(f"[panel d] operon-body mean depth (x mean): Syn1 {s1v:.3f}x, Syn3A {s3v:.3f}x, FC {s3v / s1v:.3f}")
    print(f"[panel d] mRNA-pool share: Syn1 {sh1:.2f}% -> Syn3A {sh3:.2f}%; "
          f"-10 box = {m10['minus10_6mer']} ({m10['motif_tier']})  (-> {depth_txt})")
    return out


def write_r6_stats(out_name="R6_stats.txt"):
    """Record all the key R6 numbers (panels a-e + the L6.1-L6.4 prose) into one txt."""
    df  = pd.read_csv(os.path.join(CRP, "syn1_vs_syn3a_RNA_protein.tsv"), sep="\t")
    sec = pd.read_csv(os.path.join(CRP, "TPM_change_by_secondary.tsv"), sep="\t")
    cx  = pd.read_csv(os.path.join(CRP, "macromolecule_complex_abundance.tsv"), sep="\t")
    occ = open(os.path.join(CRP, "deleted_gene_occupancy.txt")).read()
    L = []

    def fc(loc):
        r = df[df["locus_syn1"] == f"MMSYN1_{loc}"]
        if not len(r):
            return ("?", float("nan"), float("nan"))
        r = r.iloc[0]
        return (str(r["gene_name"]), float(r["TPM_fold_change"]), float(r["iPM_fold_change"]))

    L.append("R6 - Transcriptome & proteome reallocation in Syn3A : key numbers")
    L.append("=" * 64)

    L.append("\n[L6.1] Deleted-gene occupancy (raw Syn1 share)")
    for key in ("deleted loci ", "mRNA     :", "pseudo   :", "ncRNA    :", "tRNA     :",
                "mRNA pool (coding only)", "total proteome"):
        for ln in occ.splitlines():
            if key in ln:
                L.append("  " + ln.strip()); break

    mrna = df[df["rna_type"] == "mRNA"]
    both = mrna[(mrna["relTPM_syn1"] > 0) & (mrna["relTPM_syn3a"] > 0)]
    import numpy as _np
    _r, _ = stats.pearsonr(_np.log10(both["relTPM_syn1"]), _np.log10(both["relTPM_syn3a"]))
    _below = float((both["relTPM_syn3a"] < both["relTPM_syn1"]).mean() * 100)
    L.append("\n[L6.1b / panel b] Retained-pool conservation vs downshift")
    L.append(f"  relTPM Syn1 vs Syn3A  Pearson(log10) r = {_r:.3f}  (n = {len(both)} mRNA)")
    L.append(f"  genes below the y=x diagonal (relTPM_syn3A < relTPM_syn1): {_below:.1f}%  "
             "(most of the retained pool shifted down to fund the r-protein rise)")

    L.append("\n[L6.2 / panels a,b] Function-category mRNA-pool share change "
             "(retained pool, deletion-corrected; FC=Syn3A/Syn1)")
    for _, r in sec.sort_values("pool_share_change", ascending=False).iterrows():
        L.append(f"  {str(r['category'])[:26]:26s} {r['syn1_pool_share_pct']:6.2f}% -> "
                 f"{r['syn3a_pool_share_pct']:6.2f}%  (d {r['pool_share_change']:+.2f} pts)  "
                 f"medFC {r['median_TPM_FC_corr']:.2f}  p={r['mwu_p_vs_rest']:.2e}")

    L.append("\n[L6.3 / panel c] Macromolecular complexes (limiting-subunit, FC=Syn3A/Syn1)")
    for _, r in cx.iterrows():
        L.append(f"  {str(r['complex']):16s} {str(r['formula']):26s} "
                 f"TPM_FC {r['TPM_fold_change']:.2f}  iPM_FC {r['iPM_fold_change']:.2f}")
    L.append("  cell cycle: Syn3A ~105 min vs Syn1 ~60 min")

    L.append("\n[panel c] Transcript / protein FC (Syn3A/Syn1) by group")
    groups = [("PTS (Carbohydrate transport)", ["0233", "0694", "0234", "0779"]),
              ("Glycolysis (pathway order)", ["0445", "0220", "0131", "0727", "0607",
                                              "0451", "0606", "0729", "0213", "0221"]),
              ("Pyruvate metabolism", ["0227", "0229", "0230", "0475"])]
    for gname, locs in groups:
        L.append(f"  {gname}:")
        for loc in locs:
            n, t, i = fc(loc)
            L.append(f"     {n:6s}/{loc}  TPM {t:.3f}  iPM {i:.3f}")

    L.append("\n[panels d,e] Giant rPtn operon OP_00341 (MMSYN1_0652-0672, 21 genes, ~11 kb, minus strand)")
    loci = [f"MMSYN1_{n:04d}" for n in range(652, 673)]
    g = df[df["locus_syn1"].isin(loci)]
    t1, t3 = df["relTPM_syn1"].sum(), df["relTPM_syn3a"].sum()
    s1, s3 = g["relTPM_syn1"].sum(), g["relTPM_syn3a"].sum()
    L.append(f"  coding mRNA-pool share: Syn1 {100*s1/t1:.1f}%  ->  Syn3A {100*s3/t3:.1f}%  "
             f"(share FC {(s3/t3)/(s1/t1):.2f}; per-gene relTPM FC {s3/s1:.2f})")
    L.append("  single polycistron, NO internal terminator; full-length ~11kb reads rare (1-2);")
    L.append("  depth = 5'-polarity gradient (~90k at 5') with a sharp internal endonucleolytic step (tx~2100)")

    L.append("\n[L6.2 caveat] r-proteins that buck the up-trend (TPM down, protein held/up)")
    for loc, note in [("0082", "new_promoter_fusion (weak fused promoter)"),
                      ("0294", "new_promoter_fusion (weak fused promoter)"),
                      ("0526", "UNAFFECTED (operon intact, not decapitated)"),
                      ("0137", "UNAFFECTED (operon intact, not decapitated)"),
                      ("0482", "UNAFFECTED (operon intact, not decapitated)")]:
        n, t, i = fc(loc)
        L.append(f"  {n:6s}/{loc}  TPM {t:.3f}  iPM {i:.3f}  - {note}")

    L.append("\n[panel e] tRNA operon -> rPtn junction (Syn3A); deletion changed neighbour, not regulation")
    L.append("  new upstream neighbour: tRNAs MMSYN1_0678-0681 = Thr/Val/Glu/Asn (TGT/TAC/TTC/GTT), co-directional (minus)")
    L.append("  TSS(806176) -> nearest deletion (DEL_074 @806355): 179 bp  -> rPtn promoter retained/intact")
    L.append("  Syn1 TSS->tRNA-3' = 7193 bp ; deleted (DEL_074 5509 + DEL_075 912) = 6421 bp ; Syn3A = 772 bp")
    L.append("  Syn3A coords: rpsJ/0672 419409-419717 ; tRNA cluster 420560-420889 (minus)")
    L.append("  Co-expression test (ONT + Illumina, the 06/07 method) for rpsJ/0672 <-> tRNA cluster:")
    L.append("    ONT      : 0 / 3084 minus-strand reads span (enclose both); 0 bridge (>=50 bp on both genes)")
    L.append("    Illumina : true inter-operon middle (419784-420350) mean depth 27 = 1.2% of flanking "
             "(vs rpsJ 20272, tRNA 2300)")
    L.append("    -> BOTH platforms silent across the middle ; VERDICT = SPLIT (not co-transcribed)")
    L.append("  => the deletion changed the operon's neighbour but NOT its regulation; the two stay independent")

    out = os.path.join(OUTDIR, out_name)
    open(out, "w").write("\n".join(L) + "\n")
    print(f"[stats] wrote {out}")
    return out


PANELS = {"b": panel_b, "c": panel_c, "d": panel_d, "stats": write_r6_stats}

if __name__ == "__main__":
    import sys
    for k in (sys.argv[1:] or list(PANELS)):
        PANELS[k]() if k in PANELS else print(f"[skip] {k} not implemented")
