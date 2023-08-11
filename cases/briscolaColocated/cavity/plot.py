import numpy as np
from matplotlib import pyplot as plt
import os, sys

sys.path.append('../../../scripts')
import figStyle as fs

fs.prep(plt)

x = np.loadtxt('postProcessing/x/30/sample.txt')
y = np.loadtxt('postProcessing/y/30/sample.txt')

xg = np.loadtxt('../../../data/Ghia/V.txt')
yg = np.loadtxt('../../../data/Ghia/U.txt')

fig = plt.figure()

plt.plot(x[:,0], x[:,5], label='$u_y(x)$')
plt.plot(y[:,1], y[:,4], label='$u_x(y)$')

plt.plot(xg[:,0], xg[:,3], 'ok', ms=5, label='Ghia et al.')
plt.plot(yg[:,0], yg[:,3], 'ok', ms=5)

plt.xlabel('$x$, $y$')
plt.ylabel('$u_x$, $u_y$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('plot.pdf')
