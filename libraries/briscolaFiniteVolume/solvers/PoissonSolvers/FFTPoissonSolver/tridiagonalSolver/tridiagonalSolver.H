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

template<class SType>
class FFTPoissonSolver;

namespace FFT
{

template<class SType>
class tridiagonalSolver
{
    // FFT solver
    FFTPoissonSolver<SType>& solver_;

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
    const PtrList<PartialList<scalar>>& cellSizes_;

    // Diagonals
    scalarBlock D_;
    scalarList DU_;
    scalarList DL_;

    // Compute FFT eigenvalues
    void computeEigenvalues();

    // Compute tridiagonal solver diagonals
    void computeDiagonals();

public:

    // Construct from FFT solver
    tridiagonalSolver
    (
        FFTPoissonSolver<SType>& solver,
        const labelVector& BC
    );

    // Destructor
    virtual ~tridiagonalSolver();

    // Solve tridiagonal system
    void solve
    (
        scalarBlock& p,
        const scalarBlock& f,
        const bool ddt,
        const label start = 0
    );

    // Solve cyclic tridiagonal system
    void solveCyclic
    (
        scalarBlock& p,
        const scalarBlock& f,
        const bool ddt
    );

    // Solve x, y or z-direction of a scalarBlock
    void solve(scalarBlock& xyzPencil, const bool ddt);

};

}

}

}

}

#endif
