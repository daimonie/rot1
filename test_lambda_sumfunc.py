import numpy as np

x = np.arange(-1,1, 0.02)
y = np.arange(-1,1, 0.02)
k = np.arange(0,12*np.pi, 0.01*np.pi)

# define f_k = (x + y)^k
fk = lambda xx, yy, kk: ( np.cos(kk*2*np.pi/xx)**2+np.sin(kk*2*np.pi/yy)**2)**0.5
#Then,

X, Y, K = np.meshgrid(x, y, k)
# sum over k after evaluating f_k
f = fk(X, Y, K).sum(axis=-1)/(k.shape[0])

f.shape
# (100, 100)
#Finally,

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X[...,0], Y[...,0], f)
plt.show()
 