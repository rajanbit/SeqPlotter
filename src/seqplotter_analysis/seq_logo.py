import numpy as np
from collections import Counter
from math import log2
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties

# Function to calculate bit score
def _calculate_bit_score(data, bases):
	nt_scores = []
	seq_list = [list(rec.seq) for rec in data]
	seq_array = np.array(seq_list).T
	for idx, nts in enumerate(seq_array):
		nts = Counter(nts)
		nts = {nt:round(count/sum(nts.values()), 2) for nt, count in nts.items()}
		I = log2(4.0)+sum({nt:p*log2(p) for nt, p in nts.items()}.values())
		bit_score = {nt:round(p*I, 2) for nt, p in nts.items()}
		bit_score = {nt:score for nt, score in bit_score.items() if nt in bases}
		bit_score = sorted(bit_score.items(), key=lambda x: x[1])
		nt_scores.append(bit_score)
	return nt_scores

# Function to generate base info
def _base_info(base, color, property, x, y, yscale=1, ax=None):
	globalscale = 1.35
	text = property[base]
	t = mpl.transforms.Affine2D().scale(1*globalscale, yscale*globalscale) + \
		mpl.transforms.Affine2D().translate(x,y) + ax.transData
	p = PathPatch(text, lw=0, fc=color[base],  transform=t)
	if ax != None:
		ax.add_artist(p)
	return p

# Base properties for plotting seq logo | DNA
nts = {nt : TextPath((-0.30, 0), nt, size=1, prop=FontProperties(weight="bold")) for nt in "ATGC"}
nt_clr = {'G': 'orange', 'A': 'limegreen', 'C': 'blue', 'T': 'red'}

# Base properties for plotting seq logo | PROTEIN
aas = {aa : TextPath((-0.30, 0), aa, size=1, prop=FontProperties(weight="bold")) for aa in "ACDEFGHIKLMNPQRSTWYVOU"}
aa_clr = {'C': 'limegreen', 'G': 'limegreen','S': 'limegreen','T': 'limegreen','Y': 'limegreen',
	'Q': 'purple','N': 'purple', 'K': 'blue','R': 'blue','H': 'blue', 'D': 'red','E': 'red',
	'A':'black','V':'black','L':'black','I':'black','P':'black','W':'black','F':'black',
	'M':'black','O':'yellow','U':'yellow'}

