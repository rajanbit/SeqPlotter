from collections import Counter
from collections import defaultdict
import numpy as np
import pandas as pd

# Function for converting sequence object to feature matrix
def seq2feature(seq_data, k_size=3, min_seq=1):
	temp_data = {}
	var = []
	for record in seq_data:
		seq, seqid = record.seq, record.seqid
		kmers = list(map(lambda i: seq[i:i+k_size], range(len(seq) - k_size + 1)))
		temp_data[seqid] = Counter(kmers)
		var += kmers
	var = Counter(var)
	var = Counter(i for i in var.elements() if var[i] > min_seq)
	var = list(var.keys())
	mat = np.array([[counter_obj[i] if i in counter_obj else 0 for i in var] for counter_obj in temp_data.values()])
	return {"kmers":var, "sample":list(temp_data.keys()), "matrix":mat}
	#return pd.DataFrame(data=data["matrix"], index=data["sample"], columns=data["kmers"])

# Function for performing PCA on feature matrix
def PCA(seq_matrix, n_comp=2):
	x = (np.array(seq_matrix["matrix"])).astype('float64')
	x -= np.mean(x, axis = 0)
	cov = np.cov(x, rowvar = False)
	evals , evecs = np.linalg.eig(cov)
	idx = np.argsort(evals)[::-1]
	evecs = evecs[:,idx]
	evals = evals[idx]
	pcs = np.dot(x, evecs)
	return {"components":[f"PC-{i+1}" for i in range(n_comp)], "sample":seq_matrix["sample"], "matrix":np.real(pcs[:,:n_comp])}
	#return pd.DataFrame(data=np.real(pcs[:,:n_comp]), index=seq_matrix["sample"], columns=[f"PC-{i+1}" for i in range(n_comp)])
