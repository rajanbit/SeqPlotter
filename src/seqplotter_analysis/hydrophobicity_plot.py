import matplotlib.pyplot as plt

# Kyte-Doolittle scale
kydo = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
	'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
	'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
	'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }

x_plot = []
y_plot = []

# Function for converting sequence to Kyte-Doolittle propensity
def _seq_to_kydo(seq):
	kydo_values = []
	for aa in seq:
		kydo_values.append(kydo[aa])
	return(kydo_values)

# Function for smoothing the data 
def _smoothing(values_list, window):
	half_window = int((window-1)/2)
	new_values = [0]*half_window+values_list+[0]*half_window
	y = [] # Smoothened Kyte-Doolittle Values
	x = [] # Amino acid positions
	for i in range(half_window,len(new_values)-half_window):
		y.append(sum(new_values[i-half_window:i+1+half_window])/window)
	for j in range(1, len(values_list)+1):
		x.append(j)
	return(x, y)

def hydrophobicity_plot(seq, window_size=9):
	x_plot, y_plot = _smoothing(_seq_to_kydo(seq), window_size)[0], _smoothing(_seq_to_kydo(seq), window_size)[1]
	plt.plot(x_plot, y_plot)
	plt.xlabel("Amino Acid Position")
	plt.ylabel("Hydrophobicity Score")
	plt.show()
