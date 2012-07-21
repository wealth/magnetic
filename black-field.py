# -*- coding: utf-8 -*-
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
import os, sys

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#draw a vector
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

stage = sys.argv[1] if len(sys.argv) > 1 else 'initial'

data = open('data/' + stage + '.txt', 'r')

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

for i in range(len(top_xs)):
    a = Arrow3D([top_xs[i],top_xs[i]+top_us[i]],
    [right_ys[i],right_ys[i]+right_vs[i]],
    [left_zs[i],left_zs[i]+left_ws[i]], mutation_scale=5,
    lw=1, arrowstyle="-|>")
    ax.add_artist(a)

ax.auto_scale_xyz([min(top_xs+top_us), max(top_xs+top_us)],
 [min(right_ys+right_vs), max(right_ys+right_vs)],
 [min(left_zs+left_ws), max(left_zs+left_ws)])
fig.savefig('images/vectorfield-' + stage + '.png')