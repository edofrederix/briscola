import numpy as np

systems = [
    'f-colocated-scalar',
    'f-colocated-vector',
    'f-staggered-scalar_0',
    'f-staggered-scalar_1',
    'f-staggered-scalar_2',
    'f-staggered-vector_0',
    'f-staggered-vector_1',
    'f-staggered-vector_2',
]

typeSizes = [1, 3, 1, 1, 1, 3, 3, 3]

def removeBrackets(fileName):

    with open(fileName, 'r') as inFile, open(fileName + '_mod','w') as outFile:

        data = inFile.read()
        data = data.replace("(", "")
        data = data.replace(")", "")

        outFile.write(data)

for i,system in enumerate(systems):

    typeSize = typeSizes[i]

    removeBrackets(system)

    data = np.loadtxt(system + '_mod')

    A = np.array(data[:,:-typeSize])
    b = np.array(data[:,-typeSize:])

    x = np.linalg.solve(A,b)

    removeBrackets(system + '_solution')

    solution = np.loadtxt(system + '_solution_mod')

    if len(solution.shape) == 1:
        solution = solution[:, np.newaxis]

    if np.amax(abs(solution-x)) > 1e-6:
        print('Direct solver test failed for', system)
