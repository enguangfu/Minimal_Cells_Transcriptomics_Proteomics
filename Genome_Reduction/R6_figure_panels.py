#!/usr/bin/env python
"""
R6_figure_panels.py  --  born-at-size panels for Results figure R6
("Minimization reallocates transcription toward translation").

Panel inventory (R6 caption):
  a  mRNA-pool composition by Secondary function  -> reuse 09's
        Compare_RNA_Protein/mRNA_pool_composition_by_secondary.pdf  (already born-at-size)
  b  Tertiary-function share-change dumbbell      -> reuse 09's
        Compare_RNA_Protein/tertiary_share_change_dumbbell.pdf      (already born-at-size)
  c  Transcript + protein fold change of RNA polymerase, the degradosome, and the
     central-carbon enzymes  [THIS SCRIPT]
  d  Predicted ATP/GTP flux  -> blocked (needs a metabolic-flux model)

Panel c reads macromolecule_complex_abundance.tsv (RNAP, degradosome limiting-subunit
estimates) + syn1_vs_syn3a_RNA_protein.tsv (enzyme fold changes); recomputes nothing.
"""
import os
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams.update({
    "font.size": 7, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
    "ytick.labelsize": 6, "legend.fontsize": 6, "pdf.fonttype": 42, "ps.fonttype": 42,
})

GR     = os.path.dirname(os.path.abspath(__file__))
CRP    = os.path.join(GR, "Compare_RNA_Protein")
OUTDIR = os.path.join(GR, "R6_panels")
os.makedirs(OUTDIR, exist_ok=True)

TPM_C, IPM_C = "#3182bd", "#e6550d"   # transcript = blue, protein = orange


# ===================================== panel c: complex + enzyme fold-change lollipop
def panel_c(out_name="R6c_complex_enzyme_FC.pdf", figsize=(7 / 3, 7 / 2)):
    """Horizontal lollipop: per entity, transcript (TPM) and protein (iPM) fold change
    (Syn3A/Syn1) on a log axis with a reference line at 1. RNA polymerase and the
    degradosome (limiting-subunit estimates) on top, then the central-carbon enzymes."""
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
    machines = [("RNA\npolymerase",) + cxrow("RNA polymerase"),
                ("Degradosome",) + cxrow("Degradosome")]
    groups = [("Gene Expression",        GIP, machines),
              ("Carbohydrate transport", MET, [loc_row(l) for l in PTS]),
              ("Glycolysis",             MET, [loc_row(l) for l in GLY]),
              ("Pyruvate metabolism",    MET, [loc_row(l) for l in PYR])]

    GAP = 0.9
    ys, labels, lab_colors, data, headers = [], [], [], [], []
    y = 0.0
    for hdr, fam, ents in groups:
        col = PRIM[fam]
        headers.append((y + 0.6, hdr, col))
        for lab, tpm, ipm in ents:
            ys.append(y); labels.append(lab); lab_colors.append(col)
            data.append((tpm, ipm)); y -= 1.0
        y -= GAP

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.axvline(1.0, color="0.45", ls="--", lw=0.7, zorder=1)
    for yy, (tpm, ipm) in zip(ys, data):
        ax.plot([tpm, ipm], [yy, yy], color="0.75", lw=0.8, zorder=2)
        ax.scatter([tpm], [yy], facecolor=TPM_C, edgecolor="none", s=11, zorder=3)
        ax.scatter([ipm], [yy], facecolor=IPM_C, edgecolor="none", s=11, zorder=3)
    # tertiary-function header, left-aligned in the gap above each group, family-coloured
    for hy, hdr, col in headers:
        ax.text(0.015, hy, hdr, transform=ax.get_yaxis_transform(), ha="left",
                va="center", fontsize=5, fontstyle="italic", fontweight="bold", color=col)

    ax.set_xscale("log")
    ax.set_xlim(0.15, 2.2)
    ax.set_xticks([0.25, 0.5, 1, 2])
    ax.set_xticklabels(["0.25", "0.5", "1", "2"])
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=5, linespacing=0.85)
    for tick, col in zip(ax.get_yticklabels(), lab_colors):  # colour labels by Primary family
        tick.set_color(col)
    ax.set_ylim(min(ys) - 0.7, max(h[0] for h in headers) + 0.6)
    ax.set_xlabel("Fold Change (Syn3A/Syn1)", fontsize=7)
    ax.tick_params(axis="both", length=2, pad=1.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    handles = [Line2D([0], [0], marker="o", color=TPM_C, ls="", ms=4, label="Transcript"),
               Line2D([0], [0], marker="o", color=IPM_C, ls="", ms=4, label="Protein")]
    ax.legend(handles=handles, fontsize=5, loc="lower right", frameon=False,
              handletextpad=0.2, labelspacing=0.2)

    out = os.path.join(OUTDIR, out_name)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[panel c] complex+enzyme FC lollipop -> {out}")
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


PANELS = {"c": panel_c, "stats": write_r6_stats}

if __name__ == "__main__":
    import sys
    for k in (sys.argv[1:] or list(PANELS)):
        PANELS[k]() if k in PANELS else print(f"[skip] {k} not implemented")
