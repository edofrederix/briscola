import numpy as np
from matplotlib import pyplot as plt
import os, sys
import re

maxError = 2

BRISCOLA = os.environ.get('BRISCOLA')

if BRISCOLA == None:
    print("BRISCOLA environment variable not set")
    exit(1)

sys.path.append(BRISCOLA + '/scripts')
import figStyle as fs

##

if len(sys.argv) < 2:
    print("No Re specified")
    exit(1)

if len(sys.argv) < 3:
    print("No log file specified")
    exit(1)

Re = int(sys.argv[1])
logFile = str(sys.argv[2])

Res = [100, 400, 1000, 3200, 5000, 7500, 10000]
i = Res.index(Re)+1

##

x = np.loadtxt('postProcessing/x/30/sample.txt')
y = np.loadtxt('postProcessing/y/30/sample.txt')

xg = np.flip(np.loadtxt(BRISCOLA + '/data/Ghia/V.txt'), axis=0)
yg = np.flip(np.loadtxt(BRISCOLA + '/data/Ghia/U.txt'), axis=0)

xv = np.interp(xg[:,0], x[:,0], x[:,5], left=0, right=0)
yu = np.interp(yg[:,0], y[:,1], y[:,4], left=0, right=1)

ex = np.round(100*1000*np.mean(np.square(xg[:,i] - xv))**0.5)/1000
ey = np.round(100*1000*np.mean(np.square(yg[:,i] - yu))**0.5)/1000

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

rf.write(str(ex) + '\n')
rf.write(str(ey) + '\n')
rf.write('failed\n' if ex > maxError else 'passed\n')
rf.write('failed\n' if ey > maxError else 'passed\n')
rf.write(str(timeStepCount) + '\n')
rf.write(str(np.round(float(nIters)/timeStepCount*1000)/1000) + '\n')

rf.close()

##

fs.prep(plt)

fig = plt.figure()

plt.plot(xg[:,0], xg[:,i], 'ok', ms=5, label='Ghia et al.')
plt.plot(yg[:,0], yg[:,i], 'ok', ms=5)

plt.plot(x[:,0], x[:,5], label='$u_y(x)$')
plt.plot(y[:,1], y[:,4], label='$u_x(y)$')

plt.xlabel('$x$, $y$')
plt.ylabel('$u_x$, $u_y$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('plot.pdf')
