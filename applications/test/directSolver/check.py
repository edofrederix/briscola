import numpy as np

systems = [
    'f-colocated-scalar_stencil_APLU',
    'f-colocated-scalar_stencil_Eigen',
    'f-colocated-scalar_symmStencil_APLU',
    'f-colocated-scalar_symmStencil_Eigen',
    'f-colocated-vector_stencil_APLU',
    'f-colocated-vector_stencil_Eigen',
    'f-colocated-vector_symmStencil_APLU',
    'f-colocated-vector_symmStencil_Eigen',
    'f-staggered-scalar_stencil_APLU_0',
    'f-staggered-scalar_stencil_APLU_1',
    'f-staggered-scalar_stencil_APLU_2',
    'f-staggered-scalar_stencil_Eigen_0',
    'f-staggered-scalar_stencil_Eigen_1',
    'f-staggered-scalar_stencil_Eigen_2',
    'f-staggered-vector_stencil_APLU_0',
    'f-staggered-vector_stencil_APLU_1',
    'f-staggered-vector_stencil_APLU_2',
    'f-staggered-vector_stencil_Eigen_0',
    'f-staggered-vector_stencil_Eigen_1',
    'f-staggered-vector_stencil_Eigen_2',
]

typeSizes = [1, 1, 1, 1, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3]

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
