from scipy.stats import norm, gaussian_kde, rv_continuous
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os, sys

plt.rc('text',usetex=True)

def create_plot(data, title, stage):
	fig = plt.figure()
	fig.suptitle(stage + '-' + title + '-histogram')
	ax = fig.add_subplot(111)

	print('Building histogram')
	resolution = abs(max(data) - min(data))
	resolution = resolution if resolution >= 1 else 1000
	counts, bins = np.histogram(data, bins=resolution, density=True)

	print('Least Squares Polyfit')
	binc = 0.5*(bins[1:]+bins[:-1])
	ks = np.polyfit(binc, counts, 4)	
	print(ks)
	line_square = np.polyval(ks, binc)

	histogram, = ax.plot(binc, counts, 'k.')

	if not (title == 'hc' or title == 'vol'):
		least, = ax.plot(binc, counts, 'r-', linewidth=2, label='Least Squares')
		ax.set_xscale('symlog')
	else:
		least, = ax.plot(binc, line_square, 'r-', linewidth=2, label='Least Squares')
	ax.set_xlabel(title)
	ax.set_ylabel('Probability Density')
	ax.grid(True)
	fig.legend((histogram, least), ('Probability Density', 'Approximation'))

	# fig.set_size_inches(16,9)
	fig.savefig('images/' + stage + '-' + title + '-histogram.eps')	

fl = open('data/' + sys.argv[1] + '.txt', 'r')

hxs = []
hys = []
hzs = []

vols = []
hcs = []

for line in fl:
	cm = line.split()
	if sys.argv[2] == 'hx':
		hxs.append(float(cm[6]))
	elif sys.argv[2] == 'hy':
		hys.append(float(cm[7]))
	elif sys.argv[2] == 'hz':
		hzs.append(float(cm[8]))
	elif sys.argv[2] == 'vol':
		vols.append(float(cm[11]))
	elif sys.argv[2] == 'hc':
		hcs.append(float(cm[10]))

if sys.argv[2] == 'hx':
	create_plot(hxs, 'hx', sys.argv[1])
elif sys.argv[2] == 'hy':
	create_plot(hys, 'hy', sys.argv[1])
elif sys.argv[2] == 'hz':
	create_plot(hzs, 'hz', sys.argv[1])
elif sys.argv[2] == 'vol':
	create_plot(vols, 'vol', sys.argv[1])
elif sys.argv[2] == 'hc':
	create_plot(hcs, 'hc', sys.argv[1])