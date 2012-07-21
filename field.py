import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os, time, sys

from enthought.mayavi import mlab

##fig1 = plt.figure(1)
##fig1.suptitle("X-Y Projection")
##fig2 = plt.figure(2)
##fig2.suptitle("X-Z Projection")
##fig3 = plt.figure(3)
##fig3.suptitle("Y-Z Projection")

##ax1 = fig1.add_subplot(111)
##ax2 = fig2.add_subplot(111)
##ax3 = fig3.add_subplot(111)

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

# Vectors projection:
##ax1.quiver(top_xs, top_ys, top_us, top_vs)
##
##ax1.set_xlabel("X Axis")
##ax1.set_ylabel("Y Axis")
##
##ax2.quiver(left_xs, left_zs, left_us, left_ws)
##
##ax2.set_xlabel("X Axis")
##ax2.set_ylabel("Z Axis")
##
##ax3.quiver(right_ys, right_zs, right_vs, right_ws)
##
##ax3.set_xlabel("Z Axis")
##ax3.set_ylabel("Y Axis")
##

# 3d vector field
mlab.options.offscreen = True
mlab.figure(1, bgcolor=(1, 1, 1), size=(1280, 720))
q = mlab.quiver3d(top_xs, right_ys, left_zs, top_us, right_vs, left_ws, scale_mode='none')
mlab.orientation_axes()
mlab.outline(color=(0,0,0), line_width=1.0)
a = mlab.axes(color=(0.1,0.1,0.1), nb_labels=3, line_width=1.0)
a.axes.axis_label_text_property.color = (0,0,0)
q.scene.isometric_view()

mlab.savefig('images/field-' + stage + '.png')