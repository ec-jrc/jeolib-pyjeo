import numpy as np
from copy import copy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pyjeo as pj

#2-Dimensional
# jim = pj.Jim(ncol = 5, nrow = 5, otype = 'GDT_UInt16', uniform = [0,2])

# plt.figure(figsize=(10,10))
# plt.imshow(jim.np(), cmap = 'gray')
# plt.xticks([], [])
# plt.yticks([], [])
# plt.savefig('ccops_label2d_input.png')
# plt.show()

#2D function
# labeled = pj.ccops.label(jim, pj.Jim(graph=4))
#2D method
# jim.ccops.label(pj.Jim(graph=4))

# plt.figure(figsize=(10,10))
# plt.imshow(labeled.np(), cmap = 'gray')
# plt.xticks([], [])
# plt.yticks([], [])
# plt.savefig('ccops_label2d_output.png')
# plt.show()


# labeled.ccops.labelsSetArea()
# max = int(labeled.stats.getStats('max')['max'])

# plt.figure(figsize=(10,10))
# palette = copy(plt.get_cmap('gray'))
# palette.set_under('white', 1.0)  # 1.0 represents not transparent
# plt.imshow(labeled.np(), cmap = palette)
# plt.colorbar()
# plt.xticks([], [])
# plt.yticks([], [])
# plt.savefig('ccops_labelsSetArea2d_output.png')
# plt.show()

#3-Dimensional
# jim = pj.Jim(ncol = 3, nrow = 3, nplane = 3, otype = 'GDT_UInt16', uniform = [0,2])

# fig = plt.figure(figsize=(10,10))
# ax = fig.gca(projection='3d')
# volume = pj.jim2np(jim)
# colors = np.repeat(volume[:, :, :, np.newaxis], 3, axis=3)
# ax.voxels(jim.np(), facecolors=colors, edgecolors='k')
# plt.savefig('ccops_label3d_input.png')
# plt.show()

# #3D function
# labeled = pj.ccops.label(jim, pj.Jim(graph=6))
# #3D method
# jim.ccops.label(pj.Jim(graph=6))


# labeled.pixops.convert('GDT_Float32')
# labeled /= 10.0
# fig = plt.figure(figsize=(10,10))
# ax = fig.gca(projection='3d')
# volume = pj.jim2np(labeled)
# colors = np.repeat(volume[:, :, :, np.newaxis], 3, axis=3)
# ax.voxels(labeled.np(), facecolors=colors, edgecolors='k')
# plt.savefig('ccops_label3d_output.png')
# plt.show()
