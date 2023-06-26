#ifndef tridiagonalSolver_H
#define tridiagonalSolver_H

#include "List.H"
#include "scalarList.H"
#include "vector.H"
#include "block.H"
#include "fvMesh.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class FFTPoissonSolver;

namespace FFT
{

class tridiagonalSolver
{
    // FFT solver
    FFTPoissonSolver& solver_;

    // FFT transformed directions
    label dir1_;
    label dir2_;

    // Global mesh dimensions
    labelVector N_;

    // Pencil dimensions
    labelVector Nd_;

    // Pencil starting indices
    labelVector Sd_;

    // Boundary conditions
    labelVector BC_;

    // Eigenvalues
    scalarList lambda1_;
    scalarList lambda2_;

    // Cell size arrays for each direction
    const PtrList<scalarList> cellSizes_;

    // Diagonals
    scalarBlock D_;
    scalarList DU_;
    scalarList DL_;

    // Compute FFT eigenvalues
    void computeEigenvalues();

    // Compute tridiagonal solver diagonals
    void computeDiagonals();

public:

    // Construct from dictionary
    tridiagonalSolver
    (
        FFTPoissonSolver& solver,
        labelVector BC
    );

    // Destructor
    virtual ~tridiagonalSolver();

    // Solve tridiagonal system
    void solve
    (
        scalarBlock& f,
        scalarBlock& p,
        label start = 0
    );

    // Solve cyclic tridiagonal system
    void solveCyclic
    (
        scalarBlock& f,
        scalarBlock& p
    );

    // Solve x, y or z-direction of a scalarBlock
    void solve(scalarBlock& xyzPencil);

};

}

}

}

}

#endif
