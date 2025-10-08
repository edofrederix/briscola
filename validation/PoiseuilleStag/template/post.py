import numpy as np
from matplotlib import pyplot as plt
import os, sys
import re

h = 1
mu = 0.01
dpdx = 3*mu

# Quite a big tolerance needed because of the penalization IBC
maxError = 10

BRISCOLA = os.environ.get('BRISCOLA')

sys.path.append(BRISCOLA + '/scripts')
import figStyle as fs

if BRISCOLA == None:
    print("BRISCOLA environment variable not set")
    exit(1)

##

if len(sys.argv) < 2:
    print("No log file specified")
    exit(1)

logFile = str(sys.argv[1])

##

data = np.loadtxt('postProcessing/y/200/sample.txt')

y = np.linspace(-1, 1, 128)
u = 0.5/mu*dpdx*(h**2 - y**2)

u_sim = np.interp(y, data[:,1], data[:,4])

e = np.round(100*1000*np.mean(np.square(u - u_sim))**0.5)/1000

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

rf.write(str(e) + '\n')
rf.write('failed\n' if e > maxError else 'passed\n')
rf.write(str(timeStepCount) + '\n')
rf.write(str(np.round(float(nIters)/timeStepCount*1000)/1000) + '\n')

rf.close()

##

fs.prep(plt)

fig = plt.figure()

plt.plot(y, u_sim, label='Briscola')
plt.plot(y, u, '--k', label='Analytical')

plt.xlabel('$y$')
plt.ylabel('$u_x$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('plot.pdf')
