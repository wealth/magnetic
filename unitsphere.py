import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os, time, sys

fig1 = plt.figure(1, figsize=(7,7))
fig2 = plt.figure(2)
##fig3 = plt.figure(3)
##fig3.suptitle("Y-Z Projection")

ax1 = fig1.add_subplot(111, projection='3d')
ax2 = fig2.add_subplot(111)
##ax3 = fig3.add_subplot(111)

stage = sys.argv[1] if len(sys.argv) > 2 else 'initial'

data = open('data/' + stage + '.txt', 'r')

dms = []

top_xs = []
top_ys = []
top_us = []
top_vs = []

left_xs = []
left_zs = []
left_us = []
left_ws = []

right_zs = []
right_ys = []
right_ws = []
right_vs = []

for line in data:
    cm = line.split()
    top_xs.append(float(cm[0]))
    top_ys.append(float(cm[1]))
    top_us.append(float(cm[3]))
    top_vs.append(float(cm[4]))

    left_xs.append(float(cm[0]))
    left_zs.append(float(cm[2]))
    left_us.append(float(cm[3]))
    left_ws.append(float(cm[5]))

    right_zs.append(float(cm[2]))
    right_ys.append(float(cm[1]))
    right_ws.append(float(cm[5]))
    right_vs.append(float(cm[4]))

    dms.append(float(cm[9]))

us = []
vs = []
for i in range(len(top_xs)):
    vector = [top_us[i], right_vs[i], left_ws[i]]
    vector = vector / abs(np.linalg.norm(vector))
    top_us[i], right_vs[i], left_ws[i] = (vector[0], vector[1], vector[2])
    us.append(0.5 - (np.arctan2(left_ws[i], top_us[i])/(2*np.pi)))
    vs.append(0.5 - 2.0 * np.arcsin(right_vs[i]) / (2 * np.pi))

ax1.scatter(top_us, right_vs, left_ws, s = 0.01)
fig1.suptitle("Magnetic Moments: " + stage)
fig1.savefig('images/magnetic-moments-unitsphere-' + stage + '.png')

ax2.scatter(us, vs, s = 0.01)
fig2.suptitle("Magnetic Moments UV Map: " + stage)
fig2.savefig('images/magnetic-moments-uvmap-' + stage + '.png')