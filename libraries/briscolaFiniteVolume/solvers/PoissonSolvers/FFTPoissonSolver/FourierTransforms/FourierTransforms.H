#ifndef FourierTransforms_H
#define FourierTransforms_H

#include "vector.H"
#include "block.H"
#include <fftw3.h>
#include "fvMesh.H"
#include "colocatedFields.H"
#include "decomposer.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class FourierTransforms
{
	// Mesh
    const fvMesh& fvMsh_;

    // Global Size
	labelVector N_;

    // Pencil decomposer
    decomposer& decomp_;

    // Boundary condition vector
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

    // X-pencil
    scalarBlock& xPencil_;

    // Y-pencil
    scalarBlock& yPencil_;

    // FFT normalization factor
    scalar normalization_;

public:

    // Construct from mesh
    FourierTransforms
    (
        const fvMesh& fvMsh,
        decomposer& d,
        scalarBlock& xPencil,
        scalarBlock& yPencil
    );

    // Destructor
    virtual ~FourierTransforms();

    void FFTBoundaryConditions();

    void FFTWplans();

    void fwdFFTx();

    void bwdFFTx();

    void fwdFFTy();

    void bwdFFTy();

    inline labelVector& BC() { return BC_; }

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
