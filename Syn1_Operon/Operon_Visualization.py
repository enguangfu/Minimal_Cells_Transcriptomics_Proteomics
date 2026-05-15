"""
Python script for operon viasulizations.
"""

# ============================================================
# Publication-style operon plotter (v1.2) — requested changes
#
# Changes implemented:
# 1) Depth TSV is known format: 3 columns (chrom, pos, depth). Parser locked to that.
# 2) Gene arrows reflect GFF3 gene strand (NOT transcript direction), so antisense genes show opposite arrows.
#    - X-axis is still transcription-direction for the OPERON (5'→3' left->right).
# 3) Gene labels use gene name if present (Name=...), otherwise locus_tag only.
# 4) draw_ends removed — promoters/terminators will be added as annotations separately.
# 5) No top axis genomic ticks; genomic coords are shown as text labels above each gene arrow.
# 6) Gene panel split into two rows: forward (+) strand on top, reverse (-) strand on bottom.
# 7) Depth panel is optional (PLOT_DEPTH = True/False) and placed at the bottom.
# 8) Bottom x-axis shows relative transcript coordinates (5'→3').
#
# Panels (top -> bottom):
# 1) Gene annotation — forward strand row (top) and reverse strand row (bottom)
# 2) RNA isoforms on operon strand (thickness ~ count)
# 3) [Optional] Sequencing depth on operon strand only
# ============================================================

from __future__ import annotations

import os
import math
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter, FixedLocator

# ----------------------------
# Your paths (assumes you defined these earlier in notebook)
# ----------------------------
# REQUIRED globals (you provided earlier):
# MOTHER_FOLDER, OPERONS_TSV, GFF3_FILE, GENOME_FOLDER

# Isoform clusters file
MOTHER_FOLDER = ".."
ISOFOMR_FOLDER = MOTHER_FOLDER + "/isoform_annotation"
ISOFORMS_TSV   = ISOFOMR_FOLDER + "/isoform_clusters_annotated.tsv"


# ----------------------------
# Config knobs
# ----------------------------
PAD_BP = 200
PAD_BP_FRAC = 0.01
DEPTH_SMOOTH_WIN = 1  # odd looks nice; set 1 to disable

MIN_ISOFORM_FRAC = 0.005
MAX_ISOFORMS_TO_PLOT =100
ISOFORM_THICKNESS_SCALE = .5
ISOFORM_THICKNESS_MODE = "sqrt"  # "sqrt" | "linear" | "log1p"

# Gene track aesthetics
GENE_LABEL_FONTSIZE  = 18

LABEL_FONTSIZE = 18
# --- NEW: Syn3A GenBank file path ---
SYN3A_GB_FILE = MOTHER_FOLDER + "/Genomes_Input/syn3a.gb"  # <-- CHANGE to your real path

def syn3a_locus_nums_from_genbank(gb_path: str) -> set[str]:
	"""
	Return locus_tag set from a Syn3A GenBank file.
	Requires: biopython installed (from Bio import SeqIO)
	"""
	from Bio import SeqIO

	tags: set[str] = set()
	for rec in SeqIO.parse(gb_path, "genbank"):
		for feat in rec.features:
			q = feat.qualifiers
			if "locus_tag" in q and len(q["locus_tag"]) > 0:
				locus_tag = str(q["locus_tag"][0])
				locus_num = locus_tag.split('_')[1]
				tags.add(locus_num)
	return tags

# Load Syn3A locus tags
syn3a_locus_nums = syn3a_locus_nums_from_genbank(SYN3A_GB_FILE)
print("Syn3A locus nums loaded:", len(syn3a_locus_nums))

# ----------------------------
# Utilities: coordinate transforms
# ----------------------------
@dataclass(frozen=True)
class OperonCoord:
	chrom: str
	strand: str
	opid: str
	start0: int
	end0: int

	def tx_of_genome_pos0(self, pos0: int) -> int:
		"""Map genomic pos0 -> operon-local transcript coordinate (0..len), left->right is operon transcription."""
		if self.strand == "+":
			return pos0 - self.start0
		else:
			return (self.end0 - pos0)

	@property
	def length(self) -> int:
		return self.end0 - self.start0


# ----------------------------
# Read operons / isoforms / genes
# ----------------------------
# operons = pd.read_csv(OPERONS_TSV, sep="\t")
# needed = {"chrom","strand",,"start0","end0"}
# missing = needed - set(operons.columns)
# if missing:
# 	raise ValueError(f"operons.final.tsv missing columns: {sorted(missing)}")

isoforms_all = pd.read_csv(ISOFORMS_TSV, sep="\t")



def read_genes_from_gff3(gff3_path: str, syn3a_locus_nums: set[str]) -> pd.DataFrame:
	"""
	Parse GFF3 (Syn1) and mark whether each gene exists in Syn3A.
	Keeps only 'gene' features.
	"""
	rows = []
	with open(gff3_path, "r") as f:
		for line in f:
			if not line or line.startswith("#"):
				continue
			parts = line.rstrip("\n").split("\t")
			if len(parts) != 9:
				continue
			chrom, source, ftype, start1, end1, score, strand, phase, attrs = parts
			if ftype != "gene":
				continue

			s1 = int(start1)
			e1 = int(end1)
			start0 = s1 - 1          # 1-based inclusive -> 0-based
			end0   = e1              # half-open end

			# parse attributes
			ad = {}
			for kv in attrs.split(";"):
				if "=" in kv:
					k, v = kv.split("=", 1)
					ad[k] = v

			locus = ad.get("locus_tag", "") or ""
			gene_name = ad.get("Name", "") or ""

			rows.append({
				"chrom": chrom,
				"start0": start0,
				"end0": end0,
				"strand": strand,
				"locus_tag": locus,
				"gene_name": gene_name,
				"in_syn3a": (locus.split('_')[1] in syn3a_locus_nums) if locus else False,
			})

	return pd.DataFrame(rows)

# Reload genes with Syn3A existence flag
GFF3_FILE = MOTHER_FOLDER + "/Genomes_Input/syn1.genes.gff3"
genes = read_genes_from_gff3(GFF3_FILE, syn3a_locus_nums)
print("Syn1 genes parsed:", len(genes))
print("Syn1 genes present in Syn3A:", int(genes["in_syn3a"].sum()))


# ----------------------------
# Depth reader
# ----------------------------
DEPTH_BEDGRAPH_FOLDER = MOTHER_FOLDER + "/PacBio_Processing/depth_bedgraph"
DEPTH_BEDGRAPH_FILES = {
	"plus":  DEPTH_BEDGRAPH_FOLDER + "/syn1.PacBio.FLNC.HQ.plus.bedGraph",
	"minus": DEPTH_BEDGRAPH_FOLDER + "/syn1.PacBio.FLNC.HQ.minus.bedGraph",
}

def read_depth_bedfile(path: str) -> pd.DataFrame:
	"""
	Read a 4-column bedGraph (chrom, start0, end0, depth) as depth track.
	Columns are 0-based half-open intervals, depth is the 4th column.
	"""
	df = pd.read_csv(path, sep="\t", header=None, comment="#")
	if df.shape[1] < 4:
		raise ValueError(f"Depth bedGraph {path} expected 4 columns: chrom start0 end0 depth")
	df = df.iloc[:, :4]
	df.columns = ["chrom", "start0", "end0", "depth"]
	df["start0"] = df["start0"].astype(int)
	df["end0"]   = df["end0"].astype(int)
	df["depth"]  = pd.to_numeric(df["depth"], errors="coerce").fillna(0.0)
	return df


depth_plus  = read_depth_bedfile(DEPTH_BEDGRAPH_FILES["plus"])
depth_minus = read_depth_bedfile(DEPTH_BEDGRAPH_FILES["minus"])


# ----------------------------
# Subsetting helpers
# ----------------------------
def subset_intervals(df: pd.DataFrame, chrom: str, start0: int, end0: int,
					 start_col="start0", end_col="end0") -> pd.DataFrame:
	g = df[df["chrom"].astype(str) == str(chrom)]
	return g[(g[start_col] < end0) & (g[end_col] > start0)].copy()


def smooth_series(y: np.ndarray, win: int) -> np.ndarray:
	if win <= 1:
		return y
	if win % 2 == 0:
		win += 1
	k = np.ones(win, dtype=float) / win
	return np.convolve(y, k, mode="same")


def thickness_from_count(n: float) -> float:
	if ISOFORM_THICKNESS_MODE == "sqrt":
		return ISOFORM_THICKNESS_SCALE * math.sqrt(max(0.0, n))
	if ISOFORM_THICKNESS_MODE == "log1p":
		return ISOFORM_THICKNESS_SCALE * math.log1p(max(0.0, n))
	return ISOFORM_THICKNESS_SCALE * max(0.0, n)


def make_gene_label(locus_tag: str, gene_name: str) -> str:
	"""
	Requested: show locusNum/gene name, where gene name is from Name= in GFF3.
	Interpretation:
	  - If gene_name exists: use gene_name
	  - else: use locus_tag only
	If you want 'gene_name|locus_tag', change below.
	"""
	gene_name = (gene_name or "").strip()
	locus_tag = (locus_tag or "").strip()
	locus_num = locus_tag.split('_')[1]
	if gene_name:
		return locus_num + '/' + gene_name
	return locus_num

def genome_pos0_from_tx(oc: OperonCoord, x_tx: float) -> int:
	"""
	Convert operon-local tx coordinate (float) -> genomic pos0 (int).
	For + strand: genome = start0 + x
	For - strand: genome = end0 - x
	"""
	if oc.strand == "+":
		return int(round(oc.start0 + x_tx))
	else:
		return int(round(oc.end0 - x_tx))
	
def make_genome_formatter(oc: OperonCoord):
	def _fmt(x, pos=None):
		# x is tx coordinate
		if oc.strand == "+":
			g = oc.start0 + x
		else:
			g = oc.end0 - x
		return str(int(round(g)))
	return FuncFormatter(_fmt)


# ----------------------------
# Panel: Gene boxes (FancyBboxPatch, two rows)
# ----------------------------
# GENE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#                "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"]

GENE_COLORS = {"+": "#1f77b4", "-": "#ff7f0e"} # blue for + strand, red for - strand

def draw_gene_arrows(ax, oc: OperonCoord, genes_df: pd.DataFrame):
	"""
	Draw genes on a single horizontal baseline (notebook-style layout).
	Each gene is rendered as a rectangle with a triangular head pointing
	in the direction of transcription:
	  - + strand: triangle on the right end
	  - − strand: triangle on the left end
	Gene name labels sit above the rectangle. Genomic coordinates are
	shown on a secondary x-axis below the panel.
	"""
	Y_BASE = 0.5
	RECT_H = 0.28      # rectangle body height
	TRI_H  = 0.28      # triangle base height (slightly wider than rect)
	HEAD_FRAC = 0.18   # head length as fraction of gene length
	HEAD_MIN  = 20     # minimum head length in TX units (≈ bp)

	if genes_df.empty:
		ax.text(0.5, 0.7, "No genes in interval",
				transform=ax.transAxes, ha="center", va="center",
				fontsize=GENE_LABEL_FONTSIZE, fontweight="bold")

	# Single baseline running through the panel (drawn even when no genes)
	xlim_lo, xlim_hi = ax.get_xlim()
	ax.hlines(Y_BASE, xlim_lo, xlim_hi, color="black", lw=2, zorder=1)

	if not genes_df.empty:
		genes_sorted = genes_df.sort_values("start0").reset_index(drop=True)
		for _, r in genes_sorted.iterrows():
			g0 = int(r["start0"])
			g1 = int(r["end0"])
			gstrand = str(r["strand"])
			in_syn3a = bool(r.get("in_syn3a", False))

			x0 = oc.tx_of_genome_pos0(g0)
			x1 = oc.tx_of_genome_pos0(g1)
			xleft  = min(x0, x1)
			xright = max(x0, x1)
			xcen   = (xleft + xright) / 2
			width  = xright - xleft

			head_len = max(HEAD_MIN, width * HEAD_FRAC)
			head_len = min(head_len, width)

			color = GENE_COLORS[gstrand]
			points_right = (gstrand == oc.strand)

			if points_right:
				tip_x  = xright
				base_x = xright - head_len
				verts = [
					(xleft,  Y_BASE - RECT_H / 2),
					(base_x, Y_BASE - RECT_H / 2),
					(base_x, Y_BASE - TRI_H / 2),
					(tip_x,  Y_BASE),
					(base_x, Y_BASE + TRI_H / 2),
					(base_x, Y_BASE + RECT_H / 2),
					(xleft,  Y_BASE + RECT_H / 2),
				]
			else:
				tip_x  = xleft
				base_x = xleft + head_len
				verts = [
					(xright, Y_BASE - RECT_H / 2),
					(base_x, Y_BASE - RECT_H / 2),
					(base_x, Y_BASE - TRI_H / 2),
					(tip_x,  Y_BASE),
					(base_x, Y_BASE + TRI_H / 2),
					(base_x, Y_BASE + RECT_H / 2),
					(xright, Y_BASE + RECT_H / 2),
				]

			ax.add_patch(mpatches.Polygon(
				verts, closed=True,
				facecolor=color, edgecolor="black",
				lw=0, alpha=1, zorder=2,
			))

			label = make_gene_label(str(r.get("locus_tag", "")), str(r.get("gene_name", "")))
			if label:
				ax.text(xcen, Y_BASE + TRI_H / 2 + 0.06, label,
						ha="center", va="bottom",
						fontsize=GENE_LABEL_FONTSIZE, color=color, clip_on=True)
			if not in_syn3a and label:
				ax.text(xcen, Y_BASE - TRI_H / 2 - 0.06, "×",
						ha="center", va="top",
						fontsize=GENE_LABEL_FONTSIZE, color="grey", clip_on=True)

	ax.set_ylim(0, 1)
	ax.set_yticks([])
	ax.spines["left"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["bottom"].set_visible(False)
	ax.tick_params(axis="x", which="major", bottom=False, top=False, labelbottom=False)

	# ---- Secondary x-axis below: genomic absolute coordinates ----
	ax_genome = ax.twiny()
	ax_genome.set_xlim(ax.get_xlim())
	ax_genome.xaxis.set_ticks_position("bottom")
	ax_genome.xaxis.set_label_position("bottom")
	ax_genome.spines["bottom"].set_position(("outward", 8))
	ax_genome.spines["top"].set_visible(False)
	ax_genome.spines["left"].set_visible(False)
	ax_genome.spines["right"].set_visible(False)

	# Reuse the same TX tick positions already computed for the bottom axis,
	# but format them as genomic coordinates
	ax_genome.xaxis.set_major_formatter(make_genome_formatter(oc))
	ax_genome.tick_params(axis="x", which="major",
						  top=False, labeltop=False,
						  bottom=True, labelbottom=True,
						  labelsize=GENE_LABEL_FONTSIZE - 2)
	ax_genome.set_xlabel("Genomic position (bp)", fontsize=GENE_LABEL_FONTSIZE - 1, labelpad=4)

# ----------------------------
# Panel: Depth
# ----------------------------
def draw_depth(ax, oc: OperonCoord, depth_df: pd.DataFrame, plot_s0: int, plot_e0: int, strand: str):
	if depth_df.empty:
		ax.text(0.01, 0.5, "No depth", transform=ax.transAxes, va="center")
		return

	L = plot_e0 - plot_s0
	if L <= 0:
		return

	y = np.zeros(L, dtype=float)
	# depth_df here is per-base intervals; fill by overlap
	for _, r in depth_df.iterrows():
		a0 = max(plot_s0, int(r["start0"]))
		a1 = min(plot_e0, int(r["end0"]))
		if a1 <= a0:
			continue
		y[a0 - plot_s0 : a1 - plot_s0] = float(r["depth"])

	if DEPTH_SMOOTH_WIN > 1:
		y = smooth_series(y, DEPTH_SMOOTH_WIN)

	genome_positions = np.arange(plot_s0, plot_e0, dtype=int)
	x_tx = np.array([oc.tx_of_genome_pos0(int(p)) for p in genome_positions], dtype=int)

	order = np.argsort(x_tx)
	ax.plot(x_tx[order], y[order], linewidth=3, color="grey")
	strand_text = "Forward" if strand == "+" else "Reverse"
	ax.set_ylabel(f"Depth {strand_text} Strand", fontsize=LABEL_FONTSIZE)
	ax.spines["right"].set_visible(False)
	ax.spines["top"].set_visible(False)
	
	ax.tick_params(axis="y", which="major", labelsize=LABEL_FONTSIZE - 2)

# ----------------------------
# Panel: Isoforms
# ----------------------------
def _pack_intervals_min_rows(intervals: List[Tuple[int,int]]) -> List[int]:
	"""
	Given intervals [left,right] in TX coordinates, assign each to a row index so that
	no intervals overlap within a row. Greedy algorithm: sort by left, then place
	each into the first row whose current end <= left.
	Returns list of row indices aligned with input order.
	"""
	# sort by start while remembering original indices
	order = sorted(range(len(intervals)), key=lambda i: (intervals[i][0], intervals[i][1]))
	row_ends: List[int] = []   # current end for each row
	rows_out = [0] * len(intervals)

	for i in order:
		a, b = intervals[i]
		placed = False
		for r, end in enumerate(row_ends):
			if end <= a:  # no overlap (touching allowed)
				row_ends[r] = b
				rows_out[i] = r
				placed = True
				break
		if not placed:
			row_ends.append(b)
			rows_out[i] = len(row_ends) - 1

	return rows_out


def layout_isoform_tracks(iso_df: pd.DataFrame, oc: OperonCoord, plot_s0: int, plot_e0: int) -> pd.DataFrame:
	"""
	- Take top MAX_ISOFORMS_TO_PLOT by n_reads
	- Clip to plotting window and convert to TX coordinates
	- Sort by tx_left ascending (correct spatial order for both strands)
	- Global interval packing across all isoforms so non-overlapping isoforms
	  from different TSS groups can share the same y-row
	- Style (lw, alpha) scaled per TSS group relative to group's own n_reads range
	Adds columns: tx_left, tx_right, group_id, y, lw, alpha
	"""
	if iso_df.empty:
		return iso_df

	# Step 1: keep top N by abundance
	iso = iso_df.sort_values("n_reads", ascending=False).head(MAX_ISOFORMS_TO_PLOT).copy()

	# Step 2: clip in genome space and convert to TX coordinates
	tx_lefts, tx_rights = [], []
	for _, r in iso.iterrows():
		s = int(r["start0"])
		e = int(r["end0"])
		s_clip = min(max(s, plot_s0), plot_e0)
		e_clip = min(max(e, plot_s0), plot_e0)
		x0 = oc.tx_of_genome_pos0(s_clip)
		x1 = oc.tx_of_genome_pos0(e_clip)
		tx_lefts.append(min(x0, x1))
		tx_rights.append(max(x0, x1))

	iso["tx_left"]  = tx_lefts
	iso["tx_right"] = tx_rights

	# Drop zero-length after clipping
	iso = iso[iso["tx_right"] > iso["tx_left"]].copy()
	if iso.empty:
		return iso

	# Step 3: sort by tx_left ascending — correct spatial order for both strands
	iso = iso.sort_values(["tx_left", "tx_right"], ascending=[True, True]).reset_index(drop=True)

	# Step 4: global interval packing — non-overlapping isoforms share y-rows
	intervals = list(zip(iso["tx_left"].astype(int), iso["tx_right"].astype(int)))
	rows = _pack_intervals_min_rows(intervals)
	iso["y"] = [r + 1 for r in rows]

	# Step 5: style (lw, alpha) scaled within each TSS group
	def _scale_alpha(n, nmin, nmax):
		if nmax <= nmin:
			return 0.9
		return float(0.5 + 0.45 * (n - nmin) / (nmax - nmin))

	iso["lw"]    = float("nan")
	iso["alpha"] = float("nan")

	for g_start, g in iso.groupby("tx_left", sort=True):
		nvals = g["n_reads"].astype(float).to_numpy()
		nmin, nmax = float(nvals.min()), float(nvals.max())
		iso.loc[g.index, "lw"]    = [max(0.3, thickness_from_count(float(n))) for n in nvals]
		iso.loc[g.index, "alpha"] = [_scale_alpha(float(n), nmin, nmax) for n in nvals]

	# group_id for color mapping: tx_left = 5' end in TX space → same TSS = same color
	iso["group_id"] = iso["tx_left"].astype(int)

	return iso



def draw_isoforms(ax, oc: OperonCoord, iso_df: pd.DataFrame, plot_s0: int, plot_e0: int):
	if iso_df.empty:
		ax.text(0.01, 0.5, "No isoforms", transform=ax.transAxes, va="center")
		ax.set_yticks([])
		return

	iso = layout_isoform_tracks(iso_df, oc, plot_s0, plot_e0)
	# print(iso)
	if iso.empty:
		ax.text(0.01, 0.5, "No isoforms (after filter)", transform=ax.transAxes, va="center")
		ax.set_yticks([])
		return

	# --- Color mapping: same start0 -> same color ---
	# Use a qualitative colormap and stable indexing by sorted unique starts.
	starts = sorted(iso["group_id"].unique().tolist())
	cmap = plt.get_cmap("tab10")  # good for categorical groups (no explicit colors hard-coded)
	color_map = {s: cmap(i % cmap.N) for i, s in enumerate(starts)}

	# Draw each isoform as an arrow pointing in the operon transcription direction.
	# + strand: arrow left → right; - strand: arrow right → left.
	for _, r in iso.iterrows():
		left = float(r["tx_left"])
		right = float(r["tx_right"])
		y = float(r["y"])
		lw = float(r["lw"])
		alpha = float(r["alpha"])
		col = color_map[int(r["group_id"])]

		if oc.strand == "+":
			posA, posB = (left, y), (right, y)
		else:
			posA, posB = (right, y), (left, y)

		arr = mpatches.FancyArrowPatch(
			posA, posB,
			arrowstyle='-|>',
			linewidth=lw,
			color=col,
			alpha=alpha,
			shrinkA=0, shrinkB=0,
		)
		ax.add_patch(arr)

	# Tight y-limits (dynamic spacing already handled by packing)
	y_max = float(iso["y"].max()) if len(iso) else 1.0
	ax.set_ylim(-1, y_max + 1)

	ax.set_ylabel("RNA Isoforms", fontsize=LABEL_FONTSIZE)
	ax.spines["right"].set_visible(False)
	ax.spines["top"].set_visible(False)
	ax.set_yticks([])


def filter_genes_for_plot(
	genes_df: pd.DataFrame,
	chrom: str,
	plot_s0: int,
	plot_e0: int,
	gene_subset_locus: Optional[List[str]] = None,
	gene_subset_name: Optional[List[str]] = None,
	gene_subset_span: Optional[Tuple[int, int]] = None,
) -> pd.DataFrame:
	"""
	Filter genes to:
	  (1) overlap the current plotting window [plot_s0, plot_e0)
	  (2) optionally overlap gene_subset_span
	  (3) optionally match locus_tag and/or gene_name lists
	"""
	g = genes_df[(genes_df["chrom"].astype(str) == str(chrom)) &
				 (genes_df["start0"] >= plot_s0) & (genes_df["end0"] <= plot_e0)].copy()

	if gene_subset_span is not None:
		a0, a1 = int(gene_subset_span[0]), int(gene_subset_span[1])
		g = g[(g["start0"] < a1) & (g["end0"] > a0)].copy()

	if gene_subset_locus is not None:
		keep = set(map(str, gene_subset_locus))
		g = g[g["locus_tag"].astype(str).isin(keep)].copy()

	if gene_subset_name is not None:
		keep = set(map(str, gene_subset_name))
		g = g[g["gene_name"].astype(str).isin(keep)].copy()

	return g

def get_xticklabels(left: int, right: int, strand: str):
	span = abs(right - left)

	if span <= 2000:
		step = 200
	elif span <= 5000:
		step = 500
	elif span <= 15000:
		step = 1000
	else:
		step = 2000

	small = min(left, right)
	large = max(left, right)

	# Round outward to multiples of step
	start = (small // step) * step
	if start > small:
		start -= step

	end = ((large + step - 1) // step) * step

	ticks = np.arange(start, end + step, step, dtype=int)

	# Return in transcription direction
	if strand == "+":
		return ticks
	else:
		return ticks[::-1]


def plot_one_operon(
	op: pd.DataFrame,
	save_path: str,
	dpi: int = 300,
	PLOT_DEPTH: bool = True,
	isoform_reads_threshold: int = 10,
):

	chrom  = str(op["chrom"])
	strand = str(op["strand"])
	s0 = int(op["start0"])
	e0 = int(op["end0"])
	opid = str(op["operon_id"])
	oc = OperonCoord(chrom=chrom, strand=strand, opid=opid, start0=s0, end0=e0)

	# Select strand-specific depth track for the OPERON strand
	depth_df = depth_plus if strand == "+" else depth_minus

	# Default plotting window from operon final bounds
	plot_s0 = s0 - int(PAD_BP_FRAC * (e0 - s0))
	plot_e0 = e0 + int(PAD_BP_FRAC * (e0 - s0))

	# ---- genes: read locus tags directly from operon row ----
	def _parse_loci(val) -> set:
		if pd.isna(val) or str(val).strip() == "":
			return set()
		return set(str(val).split(","))

	sense_loci     = _parse_loci(op.get("sense_gene_loci", ""))
	antisense_loci = _parse_loci(op.get("antisense_gene_loci", ""))
	all_loci       = sense_loci | antisense_loci

	if not genes.empty and all_loci:
		genes_sub = genes[genes["locus_tag"].astype(str).isin(all_loci)].copy()
	else:
		genes_sub = pd.DataFrame()

	# ---- subset depth to (possibly updated) plot window ----
	depth_sub = subset_intervals(depth_df, chrom, plot_s0, plot_e0, start_col="start0", end_col="end0")

	# ---- isoforms: subset by chrom, strand, and overlap with operon span ----
	for col in ["isoform_id", "chrom", "strand", "start0", "end0", "n_reads"]:
		if col not in isoforms_all.columns:
			raise ValueError(f"ISOFORMS_TSV missing required column: {col}")
	iso_sub = isoforms_all[
		(isoforms_all["chrom"].astype(str) == chrom) &
		(isoforms_all["strand"].astype(str) == strand) &
		(isoforms_all["start0"] < e0) &
		(isoforms_all["end0"]   > s0) &
		(isoforms_all["n_reads"] >= isoform_reads_threshold)
	].copy()

	# ---- Build figure ----
	# Panels: genes (two-strand rows) | isoforms | [optional] depth
	n_panels = 3 if PLOT_DEPTH else 2
	if PLOT_DEPTH:
		height_ratios = [1.8, 2.5, 1.5]
	else:
		height_ratios = [1.8, 2.5]

	fig = plt.figure(figsize=(24, 4 * n_panels))
	gs = fig.add_gridspec(n_panels, 1, height_ratios=height_ratios, hspace=0.3)

	ax_genes    = fig.add_subplot(gs[0, 0])              # genes (fwd top / rev bottom)
	ax_isoforms = fig.add_subplot(gs[1, 0], sharex=ax_genes)  # isoforms
	ax_depth    = fig.add_subplot(gs[2, 0], sharex=ax_genes) if PLOT_DEPTH else None

	# ---- Compute x-limits (always increasing in TX space) ----
	x_left  = oc.tx_of_genome_pos0(plot_s0 if strand == "+" else plot_e0)
	x_right = oc.tx_of_genome_pos0(plot_e0 if strand == "+" else plot_s0)
	left, right = min(x_left, x_right), max(x_left, x_right)
	ax_genes.set_xlim(left, right)

	ticks_tx = get_xticklabels(left, right, strand="+")  # tx tick positions (ascending)

	# ---- Title on genes panel ----
	strand_text = "forward (+)" if strand == "+" else "reverse (-)"
	title = f"Operon {op['operon_id']} | {chrom} {strand_text} | genome [{s0}, {e0})"
	ax_genes.set_title(title, fontsize=18, loc="center")

	# ---- Panel 1: genes — forward (+) strand on top row, reverse (-) on bottom row ----
	draw_gene_arrows(ax_genes, oc, genes_sub)

	# ---- Panel 2: isoforms ----
	draw_isoforms(ax_isoforms, oc, iso_sub, plot_s0, plot_e0)

	# ---- Panel 3 (optional): depth ----
	if PLOT_DEPTH and ax_depth is not None:
		draw_depth(ax_depth, oc, depth_sub, plot_s0, plot_e0, strand)

	# ---- x-axis: hide tick labels on all panels except the bottom one ----
	bottom_ax = ax_depth if PLOT_DEPTH else ax_isoforms
	hidden_axes = [ax_genes, ax_isoforms] if PLOT_DEPTH else [ax_genes]
	for ax in hidden_axes:
		plt.setp(ax.get_xticklabels(), visible=False)

	bottom_ax.xaxis.set_major_locator(FixedLocator(ticks_tx))
	bottom_ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos=None: str(int(round(x)))))
	bottom_ax.set_xlabel("Transcript coordinate (nt)", fontsize=18)
	bottom_ax.tick_params(axis="x", which="major",
						  top=False, labeltop=False,
						  bottom=True, labelbottom=True,
						  labelsize=GENE_LABEL_FONTSIZE - 2)
	fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
	plt.close(fig)



