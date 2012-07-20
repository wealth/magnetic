from scipy.stats import norm, gaussian_kde, rv_continuous
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os, sys

def create_plot(data, title, stage):
	fig = plt.figure()
	fig.suptitle(stage + '-' + title + '-histogram')
	ax = fig.add_subplot(111)

	resolution = abs(max(data) - min(data)) / 0.1
	resolution = resolution if resolution >= 1 else 1000
	counts, bins = np.histogram(data, bins=resolution, density=True)	

	binc = 0.5*(bins[1:]+bins[:-1])
	ks = np.polyfit(binc, counts, 4)	
	line_square = np.polyval(ks, binc)
	least, = ax.plot(binc, line_square, 'r-', linewidth=2, label='Least Squares')

	density = gaussian_kde(data)
	density.covariance_factor = lambda : .25
	density._compute_covariance()
	gaussian = ax.fill_between(binc, density(binc), color=(0.1,0.1,0.1), linewidth=2)
	
	if not (title == 'hc' or title == 'dm'):
		ax.set_xscale('symlog')
	ax.set_xlabel(title)
	ax.set_ylabel('Probability Density')
	ax.grid(True)
	fig.legend((gaussian, least), ('Gaussian KDE', 'Least Squares'))

	# fig.set_size_inches(16,9)
	fig.savefig('images/' + stage + '-' + title + '-histogram.png')	

fl = open('data/' + sys.argv[1] + '.txt', 'r')

hxs = []
hys = []
hzs = []

dms = []
hcs = []

for line in fl:
	cm = line.split()
	if sys.argv[2] == 'hx':
		hxs.append(float(cm[6]))
	elif sys.argv[2] == 'hy':
		hys.append(float(cm[7]))
	elif sys.argv[2] == 'hz':
		hzs.append(float(cm[8]))
	elif sys.argv[2] == 'dm':
		dms.append(float(cm[9]))
	elif sys.argv[2] == 'hc':
		hcs.append(float(cm[10]))

if sys.argv[2] == 'hx':
	create_plot(hxs, 'hx', sys.argv[1])
elif sys.argv[2] == 'hy':
	create_plot(hys, 'hy', sys.argv[1])
elif sys.argv[2] == 'hz':
	create_plot(hzs, 'hz', sys.argv[1])
elif sys.argv[2] == 'dm':
	create_plot(dms, 'dm', sys.argv[1])
elif sys.argv[2] == 'hc':
	create_plot(hcs, 'hc', sys.argv[1])