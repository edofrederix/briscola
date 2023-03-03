#ifndef FFTPoissonSolver_H
#define FFTPoissonSolver_H

#include "PoissonSolver.H"

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
    // Function to check if the mesh qualifies for FFTs. Will give error if not
    // all criteria are met.
    void checkMesh() const;

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
