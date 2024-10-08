import numpy as np
from matplotlib import pyplot as plt
import os, sys

sys.path.append('../../../scripts')
import figStyle as fs

fs.prep(plt)

data = np.loadtxt('data.txt')

y = np.loadtxt('../../../data/Hysing/Y.txt')
u = np.loadtxt('../../../data/Hysing/U.txt')

fig = plt.figure('y')
plt.plot(data[:,0], data[:,1], label='Briscola')
plt.plot(y[:,0], y[:,1], 'ok', ms=5, mew=0.5, mec='w', label='Hysing et al.')

fig = plt.figure('u')
plt.plot(data[:,0], data[:,2], label='Briscola')
plt.plot(u[:,0], u[:,1], 'ok', ms=5, mew=0.5, mec='w', label='Hysing et al.')

#

fig = plt.figure('y')

plt.xlabel('$t$')
plt.ylabel('$y$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('y.pdf')

#

fig = plt.figure('u')

plt.xlabel('$t$')
plt.ylabel('$u$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('u.pdf')
