import numpy as np
from matplotlib import pyplot as plt
import os, sys

sys.path.append('../../../scripts')
import figStyle as fs

fs.prep(plt)

times = ['0.2', '0.4', '0.6', '0.8', '1']

fig = plt.figure()

for time in times:

    line = np.loadtxt('postProcessing/line/'+time+'/sample.txt')

    plt.plot(line[:,2], line[:,3], label='$t='+time+'$')

plt.xlabel('$z$')
plt.ylabel('$T$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('plot.pdf')
