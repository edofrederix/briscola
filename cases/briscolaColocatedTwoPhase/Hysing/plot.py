import numpy as np
from matplotlib import pyplot as plt
import os, sys

sys.path.append('../../../scripts')
import figStyle as fs


data = np.loadtxt('hysingResults.txt')

fig = plt.figure()
fs.prep(plt)
plt.rc('figure', figsize=(5, 4))

plt.plot(data[:,0], data[:,1], color = "blue")

plt.xlabel('$t$')
plt.ylabel('$y$')
plt.grid()

plt.savefig('hysingPos.pdf')

fig = plt.figure()
fs.prep(plt)
plt.rc('figure', figsize=(5, 4))

plt.plot(data[:,0], data[:,2], color = "blue")

plt.xlabel('$t$')
plt.ylabel('$v$')
plt.grid()

plt.savefig('hysingVel.pdf')