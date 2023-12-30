#ifndef FourierTransforms_H
#define FourierTransforms_H

#include "vector.H"
#include "block.H"
#include <fftw3.h>
#include "fvMesh.H"
#include "colocatedFields.H"
#include "pencilDecomposer.H"
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
class FourierTransforms
{
    // Solver
    FFTPoissonSolver<SType>& solver_;

    // Solved field
    colocatedScalarField& x_;

    // Global Size
	labelVector N_;

    // Boundary condition vector
    // 0: Empty - Empty
    // 1: Dirichlet - Dirichlet
    // 2: Neumann - Neumann
    // 3: Dirichlet - Neumann
    // 4: Neumann - Dirichlet
    // 5: Cyclic - Cyclic
    labelVector BC_;

    // FFTW plans
    fftw_plan fwdPlanX_;
    fftw_plan bwdPlanX_;
    fftw_plan fwdPlanY_;
    fftw_plan bwdPlanY_;
    fftw_plan fwdPlanZ_;
    fftw_plan bwdPlanZ_;

    // FFT normalization factor
    scalar normalization_;

public:

    // Construct from solver
    FourierTransforms
    (
        FFTPoissonSolver<SType>& solver,
        colocatedScalarField& x
    );

    // Destructor
    virtual ~FourierTransforms();

    // Set boundary condition vector
    void FFTBoundaryConditions();

    // Initialize FFTW plans
    void FFTWplans();

    // Forward and backward FFT's
    void fwdFFTx();
    void bwdFFTx();
    void fwdFFTy();
    void bwdFFTy();
    void fwdFFTz();
    void bwdFFTz();

    // Normalize transformed data
    void normalize(scalarBlock&) const;

    // Return boundary condition vector
    inline const labelVector& BC() const { return BC_; }

    // Return solved field reference
    inline const colocatedScalarField& x() const { return x_; }
};

}

}

}

}

#endif
