import numpy as np
from matplotlib import pyplot as plt
import os, sys
import re

maxError_y = 1.2
maxError_u = 1.2

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

data = np.loadtxt('data.txt')

y = np.loadtxt(BRISCOLA + '/data/Hysing/1/Y.txt')
u = np.loadtxt(BRISCOLA + '/data/Hysing/1/U.txt')

y_sim = np.interp(y[:,0], data[:,0], data[:,1])
u_sim = np.interp(u[:,0], data[:,0], data[:,2])

ey = np.round(100*1000*np.mean(np.square(y[:,1] - y_sim))**0.5)/1000
eu = np.round(100*1000*np.mean(np.square(u[:,1] - u_sim))**0.5)/1000

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

rf.write(f"{ey:.6g}\n")
rf.write(f"{eu:.6g}\n")

rf.write('failed\n' if np.isnan(ey) or ey > maxError_y else 'passed\n')
rf.write('failed\n' if np.isnan(eu) or eu > maxError_u else 'passed\n')

rf.write(str(timeStepCount) + '\n')
rf.write(str(np.round(float(nIters)/timeStepCount*1000)/1000) + '\n')

rf.close()

##

fs.prep(plt)

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
