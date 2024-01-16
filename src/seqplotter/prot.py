from seqplotter.sequence import Sequence
import matplotlib.pyplot as plt

class Protein(Sequence):

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
