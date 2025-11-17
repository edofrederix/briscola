import numpy as np
from matplotlib import pyplot as plt
import os, sys
import re

maxError_ke1 = 10
maxError_ke2 = 25

maxError_pe1 = 3
maxError_pe2 = 1

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

ke1 = np.loadtxt(BRISCOLA + '/data/phaseInversion/KE_1.txt')
ke2 = np.loadtxt(BRISCOLA + '/data/phaseInversion/KE_2.txt')
pe1 = np.loadtxt(BRISCOLA + '/data/phaseInversion/PE_1.txt')
pe2 = np.loadtxt(BRISCOLA + '/data/phaseInversion/PE_2.txt')

ke1_sim = np.interp(ke1[:,0], data[:,0], data[:,1])
ke2_sim = np.interp(ke2[:,0], data[:,0], data[:,2])
pe1_sim = np.interp(pe1[:,0], data[:,0], data[:,3])
pe2_sim = np.interp(pe2[:,0], data[:,0], data[:,4])

ek1 = np.round(100*1000*np.mean(np.square(ke1[:,1] - ke1_sim))**0.5)/1000
ek2 = np.round(100*1000*np.mean(np.square(ke2[:,1] - ke2_sim))**0.5)/1000
ep1 = np.round(100*1000*np.mean(np.square(pe1[:,1] - pe1_sim))**0.5)/1000
ep2 = np.round(100*1000*np.mean(np.square(pe2[:,1] - pe2_sim))**0.5)/1000

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

rf.write(f"{ek1:.6g}\n")
rf.write(f"{ek2:.6g}\n")
rf.write(f"{ep1:.6g}\n")
rf.write(f"{ep2:.6g}\n")

rf.write('failed\n' if np.isnan(ek1) or ek1 > maxError_ke1 else 'passed\n')
rf.write('failed\n' if np.isnan(ek2) or ek2 > maxError_ke2 else 'passed\n')
rf.write('failed\n' if np.isnan(ep1) or ep1 > maxError_pe1 else 'passed\n')
rf.write('failed\n' if np.isnan(ep2) or ep2 > maxError_pe2 else 'passed\n')

rf.write(str(timeStepCount) + '\n')
rf.write(str(np.round(float(nIters)/timeStepCount*1000)/1000) + '\n')

rf.close()

##

fs.prep(plt)

#

fig = plt.figure('ke1')
plt.plot(data[:,0], data[:,1], label='Briscola')
plt.plot(ke1[:,0], ke1[:,1], 'ok', ms=5, mew=0.5, mec='w', label='Saeedipour et al.')
plt.xlabel('$t*$')
plt.ylabel('$KE$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('ke1.pdf')

#

fig = plt.figure('ke2')
plt.plot(data[:,0], data[:,2], label='Briscola')
plt.plot(ke2[:,0], ke2[:,1], 'ok', ms=5, mew=0.5, mec='w', label='Saeedipour et al.')
plt.xlabel('$t*$')
plt.ylabel('$KE$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('ke2.pdf')

#

fig = plt.figure('pe1')
plt.plot(data[:,0], data[:,3], label='Briscola')
plt.plot(pe1[:,0], pe1[:,1], 'ok', ms=5, mew=0.5, mec='w', label='Saeedipour et al.')
plt.xlabel('$t*$')
plt.ylabel('$PE$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('pe1.pdf')

#

fig = plt.figure('pe2')
plt.plot(data[:,0], data[:,4], label='Briscola')
plt.plot(pe2[:,0], pe2[:,1], 'ok', ms=5, mew=0.5, mec='w', label='Saeedipour et al.')
plt.xlabel('$t*$')
plt.ylabel('$PE$')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('pe2.pdf')