import numpy as np
from matplotlib import pyplot as plt
import os, sys

sys.path.append('../../../scripts')
import figStyle as fs

fs.prep(plt)

y = np.loadtxt('postProcessing/y/20/sample.txt')

nu = 0.1
G = 1
y_an = np.linspace(-0.5,0.5,100)
u_an = [G/(2*nu)*(y+0.5)*(1-(y+0.5)) for y in y_an]

fig = plt.figure()

plt.plot(y[:,1], y[:,4], label='$u_x(y)$')
plt.plot(y_an, u_an, ':k', label='analytical')

plt.xlabel('$x$, $y$')
plt.ylabel('$u_x$, $u_y$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('plot.pdf')
