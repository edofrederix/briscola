#ifndef FFTPoissonSolver_H
#define FFTPoissonSolver_H

#include "PoissonSolver.H"
#include <fftw3.h>
#include "FFTPlanner.H"
#include "pencilDecomposer.H"
#include "tridiagonalSolver.H"
#include "FourierTransforms.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

using namespace FFT;

class FFTPoissonSolver
:
    public PoissonSolver<stencil,scalar,colocated>
{
    // Planner
    FFTPlanner FFTPlan_;

    // Pencil decomposer
    pencilDecomposer decomp_;

    // Initial data block
    scalarBlock initData_;

    // x-pencil block
    scalarBlock xPencil_;

    // y-pencil block
    scalarBlock yPencil_;

    // z-pencil block
    scalarBlock zPencil_;

    // FFT's
    autoPtr<FourierTransforms> fft_;

    // Tridiagonal solver
    autoPtr<tridiagonalSolver> tds_;

    // Funtion to check if the mesh qualifies for FFTs. Will give error if not
    // all criteria are met.
    void checkMesh(const fvMesh&);

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

    // Access functions
    inline pencilDecomposer& decomp() { return decomp_; }
    inline FFTPlanner& FFTPlan() { return FFTPlan_; }
    inline scalarBlock& xPencil() { return xPencil_; }
    inline scalarBlock& yPencil() { return yPencil_; }
    inline scalarBlock& zPencil() { return zPencil_; }
};

}

}

}

#endif
