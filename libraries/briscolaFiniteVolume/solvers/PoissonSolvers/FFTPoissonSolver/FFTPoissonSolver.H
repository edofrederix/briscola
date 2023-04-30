#ifndef FFTPoissonSolver_H
#define FFTPoissonSolver_H

#include "PoissonSolver.H"
#include <fftw3.h>
#include "planner.H"
#include "decomposer.H"
#include "tridiagonalSolver.H"
#include "FourierTransforms.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class FFTPoissonSolver
:
    public PoissonSolver<stencil,scalar,colocated>
{
    // Cell size arrays for each direction
    const PtrList<scalarList> cellSizes_;

    // Planner
    planner PoissonPlan_;

    // Initial decomposition type:
    // 0: Not decomposed
    // 1: 3D
    // 2: x-pencil
    // 3: y-pencil
    // 4: z-pencil
    // 5: xy-slab
    // 6: xz-slab
    // 7: yz-slab
    label decompType_;

    // Mesh uniformity per direction (0 = false, 1 = true)
    labelVector meshUniform_;

    // Solve direction
    label solveDir_;

    // Pencil decomposer
    decomposer decomp_;

    // Initial data block
    scalarBlock initData_;

    // x-pencil block
    scalarBlock xPencil_;

    // y-pencil block
    scalarBlock yPencil_;

    // z-pencil block
    scalarBlock zPencil_;

    // FFT's
    FourierTransforms fft_;

    // Tridiagonal solver
    tridiagonalSolver tds_;

    // Funtion to check if the mesh qualifies for FFTs. Will give error if not
    // all criteria are met.
    void checkMesh();

public:

    // Runtime type information
    TypeName("FFT");

    // Construct from dictionary and mesh
    FFTPoissonSolver
    (
        const word,
        const dictionary&,
        const fvMesh&
    );

    // Construct from mesh
    FFTPoissonSolver(const fvMesh&);

    // Solve function
    virtual void solve
    (
        colocatedScalarField&,
        const colocatedScalarField*,
        const colocatedFaceScalarField*,
        const bool = false
    );

    // Use the base class solve support functions
    using PoissonSolver<stencil,scalar,colocated>::solve;
};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#ifdef NoRepository
    #include "FFTPoissonSolver.C"
#endif

#endif
