from seqplotter.sequence import Sequence
import matplotlib.pyplot as plt
from collections import defaultdict
from seqplotter.nucl import DNA
import pandas as pd
from matplotlib.cm import get_cmap
from collections import Counter

class PROT(Sequence):

	# Constructor
	def __init__(self, seqid:str, seq:str):
		super().__init__(seqid, seq)

	# Amino acid composition for protein sequence
	def comp(self):
		prot_comp = super().comp() 
		temp = {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0,
			 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0,
			 'Y':0, 'V':0}
		for aa in temp:
			if aa in prot_comp:
				temp[aa] = prot_comp[aa]
		for k, v in temp.items():
			temp[k] = (v/sum(temp.values()))*100
		return temp
	
	# Barplot
	@staticmethod
	def barplot(comp_dict:dict, x_lab="Amino Acids", y_lab="Percentage (%)"):
		plt.bar(list(comp_dict.keys()), list(comp_dict.values()))
		plt.xlabel(x_lab)
		plt.ylabel(y_lab)
		plt.show()

	# Amino acids distribution plot | BOXPLOT
	@staticmethod
	def aa_distribution_plot(records):
		aa_dict = defaultdict(list)
		for record in records:
			comp = record.comp()
			for aa, freq in comp.items():
				aa_dict[aa].append(((comp[aa])/sum(comp.values()))*100)
		plt.boxplot(aa_dict.values(), patch_artist=True, medianprops=dict(color='red', linewidth=1.5))
		plt.xticks(range(1, 21), aa_dict.keys())
		plt.xlabel("Amino Acids")
		plt.ylabel("Percentage (%)")
		plt.show()

	# Length distribution plot | BOXPLOT
	@staticmethod
	def length_distribution_plot(records):
		DNA.length_distribution_plot(records, "Amino acids (aa)")

	# Per sequence length distribution plot | BARPLOT
	@staticmethod
	def length_plot(records):
		DNA.length_plot(records, "Length (in aa)")

	# Per sequence amino acid distribution plot | STACKED BARPLOT
	@staticmethod
	def per_sequence_comp(records):
		aa_dict = defaultdict(list)
		for record in records:
			aa_dict["SeqID"].append(record.seqid)
			comp = record.comp()
			for aa, freq in comp.items():
				aa_dict[aa].append(((comp[aa])/sum(comp.values()))*100)
		aa_df = pd.DataFrame(aa_dict)
		aa_df.plot(x='SeqID', kind='bar', stacked=True, figsize=(10,6), color=get_cmap("tab20").colors)
		plt.xlabel("Sequences")
		plt.ylabel("Percentage(%)")
		plt.legend(loc = "lower right")
		plt.show()

	# Calculate molecular weight of a protein sequence
	def mol_weight(self):
		aa_dict = {"A":71.0779, "C":103.1429,"D":115.0874,"E":129.114,"F":147.1739,"G":57.0513,
			"H":137.1393,"I":113.1576,"K":128.1723,"L":113.1576,"M":131.1961,"N":114.1026,
			"P":97.1152,"Q":128.1292,"R":156.1857,"S":87.0773,"T":101.1039,"V":99.1311,
			"W":186.2099,"Y":163.1733}
		aa_count_dict = Counter(list(self.seq))
		molwt = 18.0153
		for aa in aa_dict:
			molwt += aa_count_dict[aa]*aa_dict[aa]
		return(round(molwt/1000, 2))

