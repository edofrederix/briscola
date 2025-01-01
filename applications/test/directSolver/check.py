import os
import numpy as np

meshTypes = ['colocated', 'staggered']
dataTypes = ['scalar', 'vector']
stencilTypes = ['stencil', 'symmStencil']
solverTypes = [
    'default',
    'partialPivLU',
    'BiCGSTAB',
    'SuperLU',
    'Pardiso',
    'UmfPack'
]

def removeBrackets(fileName):

    with open(fileName, 'r') as inFile:

        data = inFile.read()
        data = data.replace("(", "")
        data = data.replace(")", "")

        open(fileName, 'w').write(data)

for meshType in meshTypes:
    for dataType in dataTypes:
        for stencilType in stencilTypes:
            for solverType in solverTypes:

                if meshType == 'staggered' and stencilType == 'symmStencil':
                    continue

                dataTypeSize = 1 if dataType == 'scalar' else 3
                numberOfDirections = 1 if meshType == 'colocated' else 3

                for i in range(0,numberOfDirections):

                    system = \
                        'f-' + \
                        meshType + '-' + \
                        dataType + '_' + \
                        stencilType + '_' + \
                        solverType + \
                        (('_' + str(i)) if numberOfDirections > 1 else '')

                    if not os.path.isfile(system):
                        continue

                    removeBrackets(system)
                    data = np.loadtxt(system)

                    A = np.array(data[:,:-dataTypeSize])
                    b = np.array(data[:,-dataTypeSize:])

                    x = np.linalg.solve(A,b)

                    removeBrackets(system + '_solution')
                    solution = np.loadtxt(system + '_solution')

                    if len(solution.shape) == 1:
                        solution = solution[:, np.newaxis]

                    if np.amax(abs(solution-x)) > 1e-3:
                        print('Direct solver test failed for', system)
