#ifndef FFTPoissonSolver_H
#define FFTPoissonSolver_H

#include "PoissonSolver.H"
#include <fftw3.h>
#include "decomposer.H"
#include "tridiagonalSolver.H"

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
    PtrList<scalarList> cellSizes_;

    // Pencil decomposer
    autoPtr<decomposer> decomp_;

    // Tridiagonal solver
    autoPtr<tridiagonalSolver> tds_;

    // Funtion to check if the mesh qualifies for FFTs. Will give error if not
    // all criteria are met.
    void checkMesh() const;

    // Function to prepare data
    void prepare();

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

}

}

}

#ifdef NoRepository
    #include "FFTPoissonSolver.C"
#endif

#endif
