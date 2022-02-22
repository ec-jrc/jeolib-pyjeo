import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pyjeo as pj

jim = pj.Jim(ncol = 3, nrow = 3, nplane = 3, otype = 'GDT_Byte', uniform = [0,2])
labeled = pj.ccops.label(jim, pj.Jim(graph=6))

fig = plt.figure(figsize=(10,10))
ax = fig.gca(projection='3d')
volume = pj.jim2np(jim)
colors = np.repeat(volume[:, :, :, np.newaxis], 3, axis=3)
ax.voxels(jim.np(), facecolors=colors, edgecolors='k')
plt.savefig('ccops_label_3d_input.png')
plt.show()

labeled.pixops.convert('GDT_Float32')
labeled /= 10.0
fig = plt.figure(figsize=(10,10))
ax = fig.gca(projection='3d')
volume = pj.jim2np(labeled)
colors = np.repeat(volume[:, :, :, np.newaxis], 3, axis=3)
ax.voxels(labeled.np(), facecolors=colors, edgecolors='k')
plt.savefig('ccops_label_3d_labeled.png')
plt.show()
