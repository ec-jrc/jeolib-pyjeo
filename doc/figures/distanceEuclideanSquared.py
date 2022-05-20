import matplotlib.pyplot as plt
import pyjeo as pj

jim = pj.Jim(nrow=30, ncol=30, otype='Byte')

jim[15-1:15+1,15-1:15+1] = 1
jim[15-1:15+1,15-1:15+1] = 1
jim[0,:] = 1
jim[-1,:] = 1
jim[:,0] = 1
jim[:,-1] = 1

distance = jim != 1
distance.ccops.distance2dEuclideanSquared()
dilated = distance < 4

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
plt.xticks([], [])
plt.yticks([], [])
ax.imshow(jim.np(), cmap = 'gray')
plt.savefig('distance_input.png')
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
plt.xticks([], [])
plt.yticks([], [])
ax.imshow(distance.np(), cmap = 'gray')
plt.savefig('distance_output.png')
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
plt.xticks([], [])
plt.yticks([], [])
ax.imshow(dilated.np(), cmap = 'gray')
plt.savefig('distance_threshold.png')
plt.show()
