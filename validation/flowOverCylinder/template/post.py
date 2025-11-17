import numpy as np
from matplotlib import pyplot as plt
import os, sys
import re

t_start = 50

freq_min = 0.1
freq_max = 0.3

##

if len(sys.argv) < 2:
    print("No log file specified")
    exit(1)

logFile = str(sys.argv[1])

##

data = np.loadtxt('data.txt')

t = data[:,0]
v = data[:,2]

i_start = np.argmax(t>t_start)
t = t[i_start:]
v = v[i_start:]

n = len(v)

window = np.hanning(n)
v_win = v*window

n_fft = n*16

fft = np.fft.rfft(v_win, n=n_fft)
freqs = np.fft.rfftfreq(n_fft, d=(t[-1] - t[0])/n)

power = np.abs(fft)**2

i_max = np.argmax(power)
freq_max_power = freqs[i_max]

##

log = open(logFile, 'r')

timeStepCount = 0
nIters = 0

p = re.compile("Solving for colocated p.*([0-9]+)$")

for line in log:

    if re.search('^Time =', line):
        timeStepCount += 1

    r = p.search(line)
    if r:
        nIters += int(r.group(1))

##

rf = open('result.txt', 'w')

rf.write(f"{freq_max_power:.6g}\n")

rf.write(
    'failed\n'
    if np.isnan(freq_max_power)
    or freq_max_power < freq_min
    or freq_max_power > freq_max
    else 'passed\n'
)

rf.write(str(timeStepCount) + '\n')
rf.write(str(np.round(float(nIters)/timeStepCount*1000)/1000) + '\n')

rf.close()
