import numpy as np
from matplotlib import pyplot as plt
import os, sys

sys.path.append('../../../scripts')
import figStyle as fs

Re = 200
U = 1
D = 1
t_start = 50
St_ref = 0.2

##

fs.prep(plt)

data = np.loadtxt('data.txt')

t = data[:,0]
v = data[:,2]

i_start = np.argmax(t>t_start)
t = t[i_start:]
v = v[i_start:]

n = len(v)

fft = np.fft.fft(v)
freqs = np.fft.fftfreq(n, d=(t[-1] - t[0])/n)

i_max = np.argmax(np.abs(fft))
print('freq_max =', freqs[i_max])

fig = plt.figure()

plt.plot(freqs[:n//2], np.abs(fft[:n//2]), label='$u_y$')

freq_ref = St_ref*U/D

plt.plot([freq_ref, freq_ref], [1e-1, 1e4], 'k--', label='Reference freq.')

##

plt.title('Shedding frequencies')

plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')

plt.xscale('log')
plt.yscale('log')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('plot.pdf')
