# -*- coding: utf-8 -*-
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
import os, sys
from matplotlib.patches import Rectangle

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

stage = sys.argv[1] if len(sys.argv) > 1 else 'initial'

data = open('data/' + stage + '.txt', 'r')

top_xs = []
left_zs = []
right_ys = []

dms = []

for line in data:
    cm = line.split()
    top_xs.append(float(cm[0]))
    left_zs.append(float(cm[2]))
    right_ys.append(float(cm[1]))

    dms.append(float(cm[9]))

dms = [dm * 1000 for dm in dms]
ax.scatter(top_xs,right_ys,left_zs,s=dms)


ax.auto_scale_xyz([min(top_xs), max(top_xs)],
 [min(right_ys), max(right_ys)],
 [min(left_zs), max(left_zs)])

sc = Rectangle((0, 0), 1, 1, fc="k")
fig.legend((sc,), ('Particles, e3',))
# fig.savefig('images/diameters-' + stage + '.png')

count = 1
while count < 360:
    ax.view_init(45, count + 44)
    fig.savefig('animate/' + stage + '-dms-scatter-animated-' + str(count).zfill(3) + '.png')
    count += 1
    print(count)