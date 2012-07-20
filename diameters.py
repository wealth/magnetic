import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os, sys

from enthought.mayavi import mlab

stage = sys.argv[1] if len(sys.argv) > 1 else 'initial'
data = open('data/' + stage + '.txt', 'r')

dms = []

top_xs = []
left_zs = []
right_ys = []

for line in data:
    cm = line.split()
    top_xs.append(float(cm[0]))
    left_zs.append(float(cm[2]))
    right_ys.append(float(cm[1]))

    dms.append(float(cm[9]))

# diameters
mlab.options.offscreen = True

mlab.figure(1, bgcolor=(1, 1, 1), size=(1280, 720))
dms = [x / 2.0 for x in dms]
q = mlab.quiver3d(top_xs, right_ys, left_zs, dms, dms, dms, mode='sphere', scale_factor=1)
mlab.orientation_axes()
mlab.outline(color=(0,0,0), line_width=1.0)
a = mlab.axes(color=(0.1,0.1,0.1), nb_labels=3, line_width=1.0)
a.axes.axis_label_text_property.color = (0,0,0)

q.scene.isometric_view()

count = 1
while count < 360:
    q.scene.camera.azimuth(1)
    mlab.draw()
    mlab.savefig('images/animate/' + stage + '-dms-animated-' + str(count).zfill(3) + '.png')
    count += 1
    print(count)
