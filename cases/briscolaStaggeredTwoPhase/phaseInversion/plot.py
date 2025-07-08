import numpy as np
from matplotlib import pyplot as plt
import os, sys

sys.path.append('../../../scripts')
import figStyle as fs

fs.prep(plt)

data = np.loadtxt('data.txt')

# Data from :

# Saeedipour, M., Vincent, S., & Estivalezes, J. L. (2021).
# Toward a fully resolved volume of fluid simulation of the
# phase inversion problem. Acta Mechanica, 232(7), 2695-2714.

# at 512^3 resolution

ke1 = np.loadtxt('../../../data/phaseInversion/KE_1.txt')
ke2 = np.loadtxt('../../../data/phaseInversion/KE_2.txt')
pe1 = np.loadtxt('../../../data/phaseInversion/PE_1.txt')
pe2 = np.loadtxt('../../../data/phaseInversion/PE_2.txt')

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
