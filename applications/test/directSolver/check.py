import os, re
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

for meshType in meshTypes:
    for dataType in dataTypes:
        for stencilType in stencilTypes:
            for solverType in solverTypes:

                if meshType == 'staggered' and stencilType == 'symmStencil':
                    continue

                dataTypeSize = 1 if dataType == 'scalar' else 3
                numberOfDirections = 1 if meshType == 'colocated' else 3

                system = \
                    'f-' + \
                    meshType + '-' + \
                    dataType + '_' + \
                    stencilType + '_' + \
                    solverType

                if not os.path.isfile(system):
                    continue

                indices = []
                i = 0

                with open(system) as file:
                    for line in file:
                        if re.search('^\d+ \d+$', line):
                            indices.append(i)
                        i = i+1

                indices.append(i)

                for i in range(0,numberOfDirections):

                    n = indices[i+1]-indices[i]-1

                    data = np.loadtxt(
                        system,
                        skiprows = indices[i] + 1,
                        max_rows = n
                    )

                    A = np.array(data[:,:-dataTypeSize])
                    b = np.array(data[:,-dataTypeSize:])

                    x = np.linalg.solve(A,b)

                    solution = np.loadtxt(
                        system + '_solution',
                        skiprows = indices[i]-i,
                        max_rows = n
                    )

                    if len(solution.shape) == 1:
                        solution = solution[:, np.newaxis]

                    if np.amax(abs(solution-x)) > 1e-3:
                        print('Direct solver test failed for', system)
