#ifndef FFTPoissonSolver_H
#define FFTPoissonSolver_H

#include <fftw3.h>

#include "PoissonSolver.H"
#include "planner.H"
#include "pencilDecomposer.H"
#include "tridiagonalSolver.H"
#include "FourierTransforms.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

/**

@brief Class that defines the FFT-based Poisson solver

This class defines a solver for Poisson equations with constant coefficients,
in which (for a three-diagonal problem) Fast Fourier (or Sine/Cosine, depending
on the boundary conditions) Transforms are applied in two spatial directions to
reduce the heptadiagonal linear system to a tridiagonal one. The tridiagonal
system is solved using the Thomas algorithm. Finally, the inverse transforms
are applied to give the real-space solution of the Poisson equation.

References:

Schumann, U., & Sweet, R. A. (1988).
Fast Fourier transforms for direct solution of Poisson's equation with
staggered boundary conditions.
Journal of Computational Physics, 75(1), 123-137.

Cifani, P., Kuerten, J. G. M., & Geurts, B. J. (2018).
Highly scalable DNS solver for turbulent bubble-laden channel flow.
Computers & Fluids, 172, 67-83.

*/

template<class SType>
class FFTPoissonSolver
:
    public PoissonSolver<SType,scalar,colocated>
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
    autoPtr<FFT::FourierTransforms<SType>> fft_;

    // Tridiagonal solver
    autoPtr<FFT::tridiagonalSolver<SType>> tds_;

    // Funtion to check if the mesh qualifies for FFTs. Will give error if not
    // all criteria are met.
    void checkMesh(const fvMesh&);

public:

    // Runtime type information
    TypeName("FFT");

    // Construct from dictionary and mesh
    FFTPoissonSolver
    (
        const dictionary&,
        const fvMesh&
    );

    // Construct from mesh
    FFTPoissonSolver(const fvMesh&);

    // Solve function
    virtual void solve
    (
        colocatedScalarField& x,
        const colocatedScalarField* bPtr,
        const colocatedFaceScalarField* lambdaPtr,
        const bool ddt = false,
        const scalar dtFrac = 1
    );

    // Use the base class solve support functions
    using PoissonSolver<SType,scalar,colocated>::solve;

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
