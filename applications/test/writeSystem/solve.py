import numpy as np

meshTypes = ['colocated', 'staggered']
numberOfDirections = [1, 3]

for i,meshType in enumerate(meshTypes):

    n = numberOfDirections[i]

    for d in range(0,n):

        appendix = '_'+str(d) if n > 1 else ''
        data = np.loadtxt('f-' + meshType + appendix)

        A = np.array(data[:,:-1])
        b = np.array(data[:,-1])

        x = np.linalg.solve(A,b)

        print(meshType + ' solution ' + str(d) + ' =')
        print(x)
