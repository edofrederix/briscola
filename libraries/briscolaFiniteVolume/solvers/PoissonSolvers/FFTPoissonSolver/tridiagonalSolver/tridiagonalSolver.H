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

class tridiagonalSolver
{
    // Mesh
    const fvMesh& fvMsh_;

    // Global mesh dimensions
    labelVector N_;

    // Z-pencil dimensions
    labelVector Nz_;

    // Z-pencil starting indices
    labelVector Sz_;

    // Boundary conditions
    labelVector BC_;

    // Eigenvalues
    scalarList lambdaX_;
    scalarList lambdaY_;

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
        const fvMesh& fvMsh,
        labelVector Nz,
        labelVector Sz,
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

    // Solve z-direction of a scalarBlock
    void solve(scalarBlock& zPencil);

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
