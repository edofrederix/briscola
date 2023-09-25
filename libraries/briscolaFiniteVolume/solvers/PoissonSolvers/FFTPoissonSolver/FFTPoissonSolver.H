#ifndef FFTPoissonSolver_H
#define FFTPoissonSolver_H

#include "PoissonSolver.H"
#include <fftw3.h>
#include "planner.H"
#include "pencilDecomposer.H"
#include "tridiagonalSolver.H"
#include "FourierTransforms.H"
#include "imSchemes.H"

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
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Planner
    FFT::planner FFTPlan_;

    // Pencil decomposer
    FFT::pencilDecomposer decomp_;

    // Initial data block
    scalarBlock initData_;

    // x-pencil block
    scalarBlock xPencil_;

    // y-pencil block
    scalarBlock yPencil_;

    // z-pencil block
    scalarBlock zPencil_;

    // FFT's
    autoPtr<FFT::FourierTransforms> fft_;

    // Tridiagonal solver
    autoPtr<FFT::tridiagonalSolver> tds_;

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
        const fvMesh&,
        immersedBoundary<scalar,colocated>*
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

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the decomposer
    inline FFT::pencilDecomposer& decomp() { return decomp_; }

    // Access to the FFT plan
    inline FFT::planner& FFTPlan() { return FFTPlan_; }

    // Access to the pencils
    inline scalarBlock& xPencil() { return xPencil_; }
    inline scalarBlock& yPencil() { return yPencil_; }
    inline scalarBlock& zPencil() { return zPencil_; }
};

}

}

}

#endif
