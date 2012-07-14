import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

from mayavi import mlab

h_resolution = 100000

##fig1 = plt.figure(1)
##fig1.suptitle("X-Y Projection")
##fig2 = plt.figure(2)
##fig2.suptitle("X-Z Projection")
##fig3 = plt.figure(3)
##fig3.suptitle("Y-Z Projection")
fig4 = plt.figure(4)
fig4.suptitle("hxs")
fig5 = plt.figure(5)
fig5.suptitle("hys")
fig6 = plt.figure(6)
fig6.suptitle("hzs")

##ax1 = fig1.add_subplot(111)
##ax2 = fig2.add_subplot(111)
##ax3 = fig3.add_subplot(111)
ax4 = fig4.add_subplot(111)
ax5 = fig5.add_subplot(111)
ax6 = fig6.add_subplot(111)
        
data = open('initial.txt', 'r')

dms = []

hxs = []
hys = []
hzs = []

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

    hxs.append(float(cm[6]))
    hys.append(float(cm[7]))
    hzs.append(float(cm[8]))

    dms.append(float(cm[9]))

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

# histograms
##ax4.hist(hxs, bins=h_resolution, normed=1, log=True)
##ax4.set_xscale('log')
##
##fig4.set_size_inches(32,9)
##fig4.set_dpi(320)
##fig4.savefig('hx.png', transparent=False, bbox_inches='tight', pad_inches=0)

##ax5.hist(hys, bins=h_resolution, normed=1, log=True)
##ax5.set_xscale('log')
##
##fig5.set_size_inches(32,9)
##fig5.set_dpi(320)
##fig5.savefig('hy.png', transparent=False, bbox_inches='tight', pad_inches=0)

##ax6.hist(hzs, bins=h_resolution, normed=1, log=True)
##ax6.set_xscale('log')
##
##fig6.set_size_inches(32,9)
##fig6.set_dpi(320)
##fig6.savefig('hz.png', transparent=False, bbox_inches='tight', pad_inches=0)

# 3d vector field
##mlab.quiver3d(top_xs, right_ys, left_zs, top_us, right_vs, left_ws)

# diameters
mlab.figure(1, bgcolor=(0, 0, 0), size=(1300, 700))
dms = [x / 2.0 for x in dms]
q = mlab.quiver3d(top_xs, right_ys, left_zs, dms, dms, dms, mode='sphere', scale_factor=1)
mlab.orientation_axes()
#mlab.savefig('diameters.png')
q.scene.x_plus_view()
mlab.savefig('dms-x-plus.png')
q.scene.x_minus_view()
mlab.savefig('dms-x-minus.png')
q.scene.y_plus_view()
mlab.savefig('dms-y-plus.png')
q.scene.y_minus_view()
mlab.savefig('dms-y-minus.png')
q.scene.z_plus_view()
mlab.savefig('dms-z-plus.png')
q.scene.z_minus_view()
mlab.savefig('dms-z-minus.png')
q.scene.isometric_view()
mlab.savefig('dms-isometric.png')
mlab.show()
