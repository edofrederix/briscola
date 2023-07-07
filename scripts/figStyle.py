import numpy as np

smallSize = 11
normalSize = 13

# Prepare plotting settings

def prep(plt):

    plt.rc('font', family='serif', size=11, serif='STIXGeneral')
    plt.rc('mathtext', fontset='stix')
    plt.rc('legend', handlelength=1.5, numpoints=1)

    plt.rc('font', size=normalSize)
    plt.rc('axes', titlesize=normalSize)
    plt.rc('axes', labelsize=normalSize)
    plt.rc('xtick', labelsize=smallSize)
    plt.rc('ytick', labelsize=smallSize)
    plt.rc('legend', fontsize=smallSize)
    plt.rc('figure', titlesize=normalSize)

    plt.rc('figure', figsize=(4, 4))

    return

# Style figure after plotting

def post(fig, l=False):

    ax = fig.gca()

    figSize = fig.get_size_inches()

    ax.set_position([0.17, 0.14, 0.8, 0.8*figSize[0]/figSize[1]])

    return

# Function to remove contour lines

def noLines(cnt):
    for c in cnt.collections:
        c.set_edgecolor("face")
