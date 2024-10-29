import numpy as np
from matplotlib import pyplot as plt
import os, sys

sys.path.append('../../../scripts')
import figStyle as fs

h = 1
mu = 0.01
dpdx = 4*mu

fs.prep(plt)

data = np.loadtxt('postProcessing/y/200/sample.txt')

fig = plt.figure()

plt.plot(data[:,1], data[:,4], label='Briscola')

y = np.linspace(-1, 1, 128)
u = 0.5/mu*dpdx*(h**2 - y**2)

plt.plot(y, u, '--k', label='Analytical')

plt.xlabel('$y$')
plt.ylabel('$u_x$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('plot.pdf')
