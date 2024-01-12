from SeqPlotter.sequence import Sequence
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

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
		plt.xlabel("GC")
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


#	def per_position_comp():
	
class RNA(Sequence):

	# Constructor
	def __init__(self, seqid:str, seq:str):
		super().__init__(seqid, seq)

	# Nucleotide composition for DNA sequence
	def comp(self):
		nt_comp = super().comp() 
		if "U" not in nt_comp and "T" in nt_comp:
			nt_comp["U"] = nt_comp.pop("T")
		temp = {"A":0, "U":0, "G":0, "C":0}
		for nt in temp:
			if nt in nt_comp:
				temp[nt] = nt_comp[nt]
		new_temp = {}
		for k, v in temp.items():
			new_temp[k] = (v/sum(temp.values()))*100
		return new_temp
	
	# GC percentage
	def gc_percent(self):
		return DNA.gc_percent(self)
	
	# GC percentage plot | BARPLOT
	@staticmethod
	def gc_plot(records):
		DNA.gc_plot(records)

	# GC distribution plot | BOXPLOT
	@staticmethod
	def gc_distribution_plot(records):
		DNA.gc_distribution_plot(records)

	# Nucleotide distribution plot | BOXPLOT
	@staticmethod
	def nt_distribution_plot(records):
		nt_dict = defaultdict(list)
		for record in records:
			comp = record.comp()
			nt_dict["A"].append(((comp["A"])/sum(comp.values()))*100)
			nt_dict["U"].append(((comp["U"])/sum(comp.values()))*100)
			nt_dict["G"].append(((comp["G"])/sum(comp.values()))*100)
			nt_dict["C"].append(((comp["C"])/sum(comp.values()))*100)
		plt.boxplot(nt_dict.values(), patch_artist=True, medianprops=dict(color='red', linewidth=1.5))
		plt.xticks([1, 2, 3, 4], nt_dict.keys())
		plt.xlabel("GC")
		plt.ylabel("Percentage (%)")
		plt.show()

