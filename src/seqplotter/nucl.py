from seqplotter.sequence import Sequence
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
from collections import Counter
from random import randint, uniform

class DNA(Sequence):

	# Constructor
	def __init__(self, seqid:str, seq:str):
		super().__init__(seqid, seq)
	
	# Nucleotide composition for DNA sequence
	def comp(self):
		nt_comp = super().comp() 
		temp = {"A":0, "T":0, "G":0, "C":0}
		for nt in temp:
			if nt in nt_comp:
				temp[nt] = nt_comp[nt]
		new_temp = {}
		for k, v in temp.items():
			new_temp[k] = (v/sum(temp.values()))*100
		return new_temp
	
	# Barplot
	@staticmethod
	def barplot(comp_dict:dict, x_lab="Nucleotides", y_lab="Percentage (%)"):
		plt.bar(list(comp_dict.keys()), list(comp_dict.values()))
		plt.xlabel(x_lab)
		plt.ylabel(y_lab)
		plt.show()

	# GC percentage
	def gc_percent(self):
		comp = self.comp()
		return f'{((comp["G"]+comp["C"])/sum(comp.values()))*100:.2f}'

	# GC percentage plot | BARPLOT
	@staticmethod
	def gc_plot(records):
		seq_dict = {}
		for record in records:
			seq_dict[record.seqid] = float(record.gc_percent())
		DNA.barplot(seq_dict, "Sequences", "GC Content (%)")
	
	# GC distribution plot | BOXPLOT
	@staticmethod
	def gc_distribution_plot(records):
		distr_lis = []
		for record in records:
			distr_lis.append(float(record.gc_percent()))
		plt.boxplot(distr_lis, patch_artist=True, medianprops=dict(color='red', linewidth=1.5))
		plt.xlabel("GC")
		plt.ylabel("Percentage (%)")
		plt.show()

	# Nucleotide distribution plot | BOXPLOT
	@staticmethod
	def nt_distribution_plot(records):
		nt_dict = defaultdict(list)
		for record in records:
			comp = record.comp()
			nt_dict["A"].append(((comp["A"])/sum(comp.values()))*100)
			nt_dict["T"].append(((comp["T"])/sum(comp.values()))*100)
			nt_dict["G"].append(((comp["G"])/sum(comp.values()))*100)
			nt_dict["C"].append(((comp["C"])/sum(comp.values()))*100)
		plt.boxplot(nt_dict.values(), patch_artist=True, medianprops=dict(color='red', linewidth=1.5))
		plt.xticks([1, 2, 3, 4], nt_dict.keys())
		plt.xlabel("Nucleotides")
		plt.ylabel("Percentage (%)")
		plt.show()
	
	# Per sequence nucleotide distribution plot | STACKED BARPLOT
	@staticmethod
	def per_sequence_comp(records):
		nt_dict = defaultdict(list)
		for record in records:
			comp = record.comp()
			nt_dict["SeqID"].append(record.seqid)
			nt_dict["A"].append(((comp["A"])/sum(comp.values()))*100)
			nt_dict["T"].append(((comp["T"])/sum(comp.values()))*100)
			nt_dict["G"].append(((comp["G"])/sum(comp.values()))*100)
			nt_dict["C"].append(((comp["C"])/sum(comp.values()))*100)		
		plt.bar(nt_dict["SeqID"], np.array(nt_dict["A"]), color='r')
		plt.bar(nt_dict["SeqID"], np.array(nt_dict["T"]), bottom=np.array(nt_dict["A"]), color='b')
		plt.bar(nt_dict["SeqID"], np.array(nt_dict["G"]), bottom=np.array(nt_dict["A"])+np.array(nt_dict["T"]), color='y')
		plt.bar(nt_dict["SeqID"], np.array(nt_dict["C"]), bottom=np.array(nt_dict["A"])+np.array(nt_dict["T"])+np.array(nt_dict["G"]), color='g')
		plt.xlabel("Sequences")
		plt.ylabel("Percentage")
		plt.legend(["A", "T", "G", "C"])
		plt.show()

	# Length distribution plot | BOXPLOT
	@staticmethod
	def length_distribution_plot(records, y_lab="Base pairs (bp)"):
		distr_lis = []
		for record in records:
			distr_lis.append(float(record.length()))
		plt.boxplot(distr_lis, patch_artist=True, medianprops=dict(color='red', linewidth=1.5))
		plt.xticks([1], ["Length"])
		plt.ylabel(y_lab)
		plt.show()

	# Per sequence length distribution plot | BARPLOT
	@staticmethod
	def length_plot(records, y_lab="Length (in bp)"):
		seq_dict = {}
		for record in records:
			seq_dict[record.seqid] = float(record.length())
		DNA.barplot(seq_dict, "Sequences", y_lab)

	# K-mer frequency distribution plot | HISTOGRAM
	def kmer_abundance_plot(self, k=3):
		kmers = [self.seq[i:i+k] for i in range(len(self.seq)-k+1)]
		kmers = Counter(kmers)
		freq = list(kmers.values())
		plt.hist(freq, bins=range(1, max(freq)+2), color="#1f77b4", edgecolor="black")
		plt.xlabel("Frequency")
		plt.ylabel("K-mers Count")
		plt.show()

	# Combined kmer frequency distribution plot | HISTOGRAM
	@staticmethod
	def combined_kmer_abundance_plot(records, k=3):
		temp = Counter()
		for record in records:
			kmers = [record.seq[i:i+k] for i in range(len(record.seq)-k+1)]
			temp += Counter(kmers)
		freq = list(temp.values())
		plt.hist(freq, bins=range(1, max(freq)+2), color="#1f77b4", edgecolor="black")
		plt.xlabel("Frequency")
		plt.ylabel("K-mers Count")
		plt.show()

	# Function to translate gene sequence (DNA) into amino acid sequence
	def translate(self):
		if len(self.seq)%3 == 0 and set(super().comp().keys()).issubset(["A", "T", "G", "C"])==True:
			codon_table = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        			'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        			'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        			'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        			'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        			'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        			'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        			'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        			'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        			'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        			'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        			'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        			'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        			'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        			'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        			'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
			codons = [self.seq[i:i+3] for i in range(0, len(self.seq)-3+1, 3)]
			return "".join([codon_table[c] for c in codons])

		elif len(self.seq)%3 != 0:
			raise IndexError("Gene sequence not in frame")

		elif set(super().comp().keys()).issubset(["A", "T", "G", "C"])==False:
			raise TypeError("Gene sequence contain ambiguous nucleotides")

	# Function to generate random DNA sequence with length=l and gc_frequency=g
	@classmethod
	def generate_random_seq(cls, l:int, gc:float, count=1):
		if count == 1:
			seq = ""
			nt = ["A", "T", "G", "C"]
			for i in range(l):
				if uniform(0.01, 1.0) < gc:
					seq += nt[randint(2, 3)]
				else:
					seq += nt[randint(0, 1)]
			cls("random_seq", seq)
		else:
			for c in range(count):
				seq = ""
				nt = ["A", "T", "G", "C"]
				for i in range(l):
					if uniform(0.01, 1.0) < gc:
						seq += nt[randint(2, 3)]
					else:
						seq += nt[randint(0, 1)]
				cls(f"random_seq{c+1}", seq)

class RNA(DNA):

	# Constructor
	def __init__(self, seqid:str, seq:str):
		super().__init__(seqid, seq)

