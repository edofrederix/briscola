#ifndef tridiagonalSolver_H
#define tridiagonalSolver_H

#include "List.H"
#include "scalarList.H"
#include "vector.H"
#include "block.H"
#include "fvMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class tridiagonalSolver
{
    // Mesh
    const fvMesh& fvMsh_;

    // Global mesh dimensions
    labelVector N_;

    // Z-pencil dimensions
    labelVector Nz_;

    // Z-pencil starting indices
    labelVector sz_;

    // Boundary conditions
    labelVector BC_;

    // Eigenvalues
    scalarList lambdaX_;
    scalarList lambdaY_;

    // Cell size arrays for each direction
    PtrList<scalarList> cellSizes_;

public:

    // Construct from dictionary
    tridiagonalSolver
    (
        const fvMesh& fvMsh,
        labelVector Nz,
        labelVector sz,
        labelVector BC
    );

    // Destructor
    virtual ~tridiagonalSolver();

    // Solve tridiagonal system
    void solve
    (
        label n,
        double DU[],
        double D[],
        double DL[],
        double f[],
        double p[]
    );

    // Solve cyclic tridiagonal system
    void solveCyclic
    (
        label n,
        double DU[],
        double D[],
        double DL[],
        double f[],
        double p[]
    );

    // Solve z-direction of a scalarBlock
    void solve(scalarBlock& zPencil);

    // Compute FFT eigenvalues
    void computeEigenvalues();

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
