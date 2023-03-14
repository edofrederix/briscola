#ifndef FFTPoissonSolver_H
#define FFTPoissonSolver_H

#include "PoissonSolver.H"
#include <fftw3.h>
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
    PtrList<scalarList> cellSizes_;

    // Pencil decomposer
    autoPtr<decomposer> decomp_;

    // Tridiagonal solver
    autoPtr<tridiagonalSolver> tds_;

    // FFT's
    autoPtr<FourierTransforms> fft_;

    // x-pencil block
    scalarBlock xPencil_;

    // y-pencil block
    scalarBlock yPencil_;

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

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#ifdef NoRepository
    #include "FFTPoissonSolver.C"
#endif

#endif
