class Sequence:
	
	records = []
	
	# Constructor
	def __init__(self, seqid:str, seq:str):
		
		# Validating sequence data
		assert seq != "", f"sequence not found for {seqid}"
		assert seqid != "", f"record id not found for {seq}"
		
		# Assigning sequence data to self variables
		self.seq = seq
		self.seqid = seqid
		
		# Adding sequence data to records
		Sequence.records.append(self)

	# Representers
	def __repr__(self):
		return f"{self.__class__.__name__}('{self.seqid}', '{self.seq}')"

	def __str__(self):
		return f"{self.seqid}:{self.seq[:70]}..."
	
	# Read FASTA/Multi-FASTA records and 
	# store it as Sequence object i.e. Sequence('SeqID', 'Sequence')
	@classmethod
	def read_fasta(cls, file_path):
		cls.records.clear()
		with open(file_path, "r") as fasta:
			seq = ""
			for line in fasta.readlines():
				if line[0] == ">" and seq == "":
					header = line[1:].strip()
				elif line[0] != ">":
					seq += line.strip()
				elif line[0] == ">" and seq != "":
					cls(seqid=header, seq=seq)
					header = line[1:].strip()
					seq = ""
			cls(seqid=header, seq=seq)
		return cls.records

	# Calculate sequence length
	def length(self):
		return len(self.seq)

	# Fetch sequence composition
	def comp(self):
		comp_dict = {}
		for s in set(self.seq):
			comp_dict[s] = self.seq.count(s)
		return comp_dict

	# Sequence count
	@staticmethod
	def count(records):
		return len(records)

	# Print head (default: first 5 lines)
	@staticmethod
	def head(nrows=5):
		out_str = "seqID\tseq\n"
		for record in Sequence.records[:nrows]:
			out_str += f"{record.seqid}:\t[{record.seq[:70]}]\n"
		print(out_str.strip())

	# Print tail (default: last 5 lines)
	@staticmethod
	def tail(nrows=5):
		out_str = "seqID\tseq\n"
		for record in Sequence.records[-(nrows):]:
			out_str += f"{record.seqid}:\t[{record.seq[:70]}]\n"
		print(out_str.strip())

	# Slice sequence
	def slice(self, start:int, end=None):
		if end == None and start <= len(self.seq):
			return self.seq[start-1:start]
		elif start <= len(self.seq) and start < end and end >=0:
			return self.seq[start-1:end]
		else:
			return("Input valid arguments (start, end) for slicing")

	# Clearing records
	@classmethod
	def clear(cls):
		cls.records.clear()

	# Finding records using seqID
	@classmethod
	def where(cls, **kwargs):
		data = {i.seqid:i.seq for i in cls.records}
		for option, value in kwargs.items():
			if option == "seqid" and value in data.keys():
				return(data[value])
			elif option == "seqid" and value not in data.keys():
				raise KeyError(f'seqID "{value}" not found')
			elif option != "seqid":
				raise ValueError(f'Invalid option: {option}')

