#ifndef tridiagonalSolver_H
#define tridiagonalSolver_H

#include "List.H"
#include "vector.H"
#include "block.H"

namespace Foam
{

namespace briscola
{

class tridiagonalSolver
{

    // System size
    label N_;

public:

    // Construct from dictionary
    tridiagonalSolver(label n);

    // Destructor
    virtual ~tridiagonalSolver();

    // Solve tridiagonal system
    void tridiagonalSolve
    (
        label n,
        double DU[],
        double D[],
        double DL[],
        double f[],
        double p[]
    );

    // Solve cyclic tridiagonal system
    void cyclicTridiagonalSolve
    (
        label n,
        double DU[],
        double D[],
        double DL[],
        double f[],
        double p[]
    );

};

}

}

#endif
