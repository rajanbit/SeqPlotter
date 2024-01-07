from SeqPlotter.sequence import Sequence
import matplotlib.pyplot as plt

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
		for k, v in temp.items():
			temp[k] = (v/sum(temp.values()))*100
		return temp
	
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
	# GC percentage plot
	@staticmethod
	def gc_plot(records):
		seq_dict = {}
		for record in records:
			seq_dict[record.seqid] = float(record.gc_percent())
		DNA.barplot(seq_dict, "Sequences", "GC Content (%)")

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
		for k, v in temp.items():
			temp[k] = (v/sum(temp.values()))*100
		return temp

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
	
	# GC percentage plot
	@staticmethod
	def gc_plot(records):
		seq_dict = {}
		for record in records:
			seq_dict[record.seqid] = float(record.gc_percent())
		RNA.barplot(seq_dict, "Sequences", "GC Content (%)")

