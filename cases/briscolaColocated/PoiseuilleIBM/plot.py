import numpy as np
from matplotlib import pyplot as plt
import os, sys

sys.path.append('../../../scripts')
import figStyle as fs

fs.prep(plt)

y = np.loadtxt('postProcessing/y/20/sample.txt')

fig = plt.figure()

plt.plot(y[:,1], y[:,4], label='$u_x(y)$')

plt.xlabel('$x$, $y$')
plt.ylabel('$u_x$, $u_y$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('plot.pdf')
