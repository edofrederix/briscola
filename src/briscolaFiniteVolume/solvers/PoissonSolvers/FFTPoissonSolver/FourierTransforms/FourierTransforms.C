#include "FourierTransforms.H"
#include "FFTPoissonSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace FFT
{

// Constructor

template<class SType>
FourierTransforms<SType>::FourierTransforms
(
    FFTPoissonSolver<SType>& solver,
    colocatedScalarField& x
)
:
    solver_(solver),
    x_(x),
    N_(solver.fvMsh().msh().template cast<rectilinearMesh>().N())
{
    FFTBoundaryConditions();
    FFTWplans();
}

// Destructor

template<class SType>
FourierTransforms<SType>::~FourierTransforms()
{
    fftw_destroy_plan(fwdPlanX_);
    fftw_destroy_plan(bwdPlanX_);
    fftw_destroy_plan(fwdPlanY_);
    fftw_destroy_plan(bwdPlanY_);
    fftw_destroy_plan(fwdPlanZ_);
    fftw_destroy_plan(bwdPlanZ_);
}

template<class SType>
void FourierTransforms<SType>::FFTBoundaryConditions()
{
    BC_ = labelVector(-1,-1,-1);

    switch (globalBoundaryConditionBaseType(x_, faceOffsets[0]))
    {
        case 1:
            BC_.x() = 0; // Empty
            break;

        case 3:
            BC_.x() = 5; // C-C
            break;

        case 4:
            if (globalBoundaryConditionBaseType(x_, faceOffsets[1]) == 4)
            {
                BC_.x() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(x_, faceOffsets[1]) == 5)
            {
                BC_.x() = 3; // D-N
            }
            break;

        case 5:
            if (globalBoundaryConditionBaseType(x_, faceOffsets[1]) == 4)
            {
                BC_.x() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(x_, faceOffsets[1]) == 5)
            {
                BC_.x() = 2; // N-N
            }
            break;

        default:
            FatalError
                << "Incorrect boundary condition." << endl
                << abort(FatalError);
            break;
    }

    switch (globalBoundaryConditionBaseType(x_, faceOffsets[2]))
    {
        case 1:
            BC_.y() = 0; // Empty
            break;

        case 3:
            BC_.y() = 5; // C-C
            break;

        case 4:
            if (globalBoundaryConditionBaseType(x_, faceOffsets[3]) == 4)
            {
                BC_.y() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(x_, faceOffsets[3]) == 5)
            {
                BC_.y() = 3; // D-N
            }
            break;

        case 5:
            if (globalBoundaryConditionBaseType(x_, faceOffsets[3]) == 4)
            {
                BC_.y() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(x_, faceOffsets[3]) == 5)
            {
                BC_.y() = 2; // N-N
            }
            break;

        default:
            FatalError
                << "Incorrect boundary condition." << endl
                << abort(FatalError);
            break;
    }

    switch (globalBoundaryConditionBaseType(x_, faceOffsets[4]))
    {
        case 1:
            BC_.z() = 0; // Empty
            break;

        case 3:
            BC_.z() = 5; // C-C
            break;

        case 4:
            if (globalBoundaryConditionBaseType(x_, faceOffsets[5]) == 4)
            {
                BC_.z() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(x_, faceOffsets[5]) == 5)
            {
                BC_.z() = 3; // D-N
            }
            break;

        case 5:
            if (globalBoundaryConditionBaseType(x_, faceOffsets[5]) == 4)
            {
                BC_.z() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(x_, faceOffsets[5]) == 5)
            {
                BC_.z() = 2; // N-N
            }
            break;

        default:
            FatalError
                << "Incorrect boundary condition." << endl
                << abort(FatalError);
            break;
    }
}

template<class SType>
void FourierTransforms<SType>::FFTWplans()
{
    fftw_r2r_kind transforms[] =
    {
        // forward transform types
        FFTW_RODFT10, FFTW_RODFT10,
        FFTW_REDFT10, FFTW_RODFT11,
        FFTW_REDFT11, FFTW_R2HC,
        // backward transform types
        FFTW_RODFT01, FFTW_RODFT01,
        FFTW_REDFT01, FFTW_RODFT11,
        FFTW_REDFT11, FFTW_HC2R
    };

    // Transform type based on boundary conditions

    fftw_r2r_kind kind_fwd_x[] = {transforms[BC_.x()]};
    fftw_r2r_kind kind_bwd_x[] = {transforms[BC_.x()+6]};

    fftw_r2r_kind kind_fwd_y[] = {transforms[BC_.y()]};
    fftw_r2r_kind kind_bwd_y[] = {transforms[BC_.y()+6]};

    fftw_r2r_kind kind_fwd_z[] = {transforms[BC_.z()]};
    fftw_r2r_kind kind_bwd_z[] = {transforms[BC_.z()+6]};

    // Pencil dimensions

    labelVector Nx = solver_.decomp().Nx()[Pstream::myProcNo()];
    labelVector Ny = solver_.decomp().Ny()[Pstream::myProcNo()];
    labelVector Nz = solver_.decomp().Nz()[Pstream::myProcNo()];

    // x-FFT plans

    int fft_rank = 1;
    int howmany = Nx.y() * Nx.z();
    int stride = 1;
    int sep = N_.x();
    const int Na_x[] = {N_.x()};

    fwdPlanX_ = fftw_plan_many_r2r
    (
        fft_rank,
        Na_x,
        howmany,
        reinterpret_cast<double*>(solver_.xPencil().begin()),
        Na_x,
        stride,
        sep,
        reinterpret_cast<double*>(solver_.xPencil().begin()),
        Na_x,
        stride,
        sep,
        kind_fwd_x,
        FFTW_MEASURE
    );

    bwdPlanX_ = fftw_plan_many_r2r
    (
        fft_rank,
        Na_x,
        howmany,
        reinterpret_cast<double*>(solver_.xPencil().begin()),
        Na_x,
        stride,
        sep,
        reinterpret_cast<double*>(solver_.xPencil().begin()),
        Na_x,
        stride,
        sep,
        kind_bwd_x,
        FFTW_MEASURE
    );

    // y-FFT plans

    fft_rank = 1;
    howmany = Ny.x() * Ny.z();
    stride = 1;
    sep = N_.y();
    const int Na_y[] = {N_.y()};

    fwdPlanY_ = fftw_plan_many_r2r
    (
        fft_rank,
        Na_y,
        howmany,
        reinterpret_cast<double*>(solver_.yPencil().begin()),
        Na_y,
        stride,
        sep,
        reinterpret_cast<double*>(solver_.yPencil().begin()),
        Na_y,
        stride,
        sep,
        kind_fwd_y,
        FFTW_MEASURE
    );

    bwdPlanY_ = fftw_plan_many_r2r
    (
        fft_rank,
        Na_y,
        howmany,
        reinterpret_cast<double*>(solver_.yPencil().begin()),
        Na_y,
        stride,
        sep,
        reinterpret_cast<double*>(solver_.yPencil().begin()),
        Na_y,
        stride,
        sep,
        kind_bwd_y,
        FFTW_MEASURE
    );

    // z-FFT plans

    fft_rank = 1;
    howmany = Nz.x() * Nz.y();
    stride = 1;
    sep = N_.z();
    const int Na_z[] = {N_.z()};

    fwdPlanZ_ = fftw_plan_many_r2r
    (
        fft_rank,
        Na_z,
        howmany,
        reinterpret_cast<double*>(solver_.zPencil().begin()),
        Na_z,
        stride,
        sep,
        reinterpret_cast<double*>(solver_.zPencil().begin()),
        Na_z,
        stride,
        sep,
        kind_fwd_z,
        FFTW_MEASURE
    );

    bwdPlanZ_ = fftw_plan_many_r2r
    (
        fft_rank,
        Na_z,
        howmany,
        reinterpret_cast<double*>(solver_.zPencil().begin()),
        Na_z,
        stride,
        sep,
        reinterpret_cast<double*>(solver_.zPencil().begin()),
        Na_z,
        stride,
        sep,
        kind_bwd_z,
        FFTW_MEASURE
    );
}

template<class SType>
void FourierTransforms<SType>::fwdFFTx()
{
    if (BC_.x()!= 0)
    {
        normalization_ = 1.0;
        fftw_execute(fwdPlanX_);
    }
}

template<class SType>
void FourierTransforms<SType>::bwdFFTx()
{
    if (BC_.x() != 0)
    {
        if (BC_.x() == 5)
        {
            normalization_ *= N_.x();
            fftw_execute(bwdPlanX_);
        }
        else
        {
            normalization_ *= 2.0 * N_.x();
            fftw_execute(bwdPlanX_);
        }
    }
}

template<class SType>
void FourierTransforms<SType>::fwdFFTy()
{
    if (BC_.y() != 0)
    {
        normalization_ = 1.0;
        fftw_execute(fwdPlanY_);
    }
}

template<class SType>
void FourierTransforms<SType>::bwdFFTy()
{
    if (BC_.y() != 0)
    {
        if (BC_.y() == 5)
        {
            normalization_ *= N_.y();
        }
        else
        {
            normalization_ *= 2.0 * N_.y();
        }
        fftw_execute(bwdPlanY_);
    }
}

template<class SType>
void FourierTransforms<SType>::fwdFFTz()
{
    if (BC_.z() != 0)
    {
        normalization_ = 1.0;
        fftw_execute(fwdPlanZ_);
    }
}

template<class SType>
void FourierTransforms<SType>::bwdFFTz()
{
    if (BC_.z() != 0)
    {
        if (BC_.z() == 5)
        {
            normalization_ *= N_.z();
        }
        else
        {
            normalization_ *= 2.0 * N_.z();
        }
        fftw_execute(bwdPlanZ_);
    }
}

template<class SType>
void FourierTransforms<SType>::normalize(scalarBlock& transformedData) const
{
    transformedData /= normalization_;
}

// Instantiate

template class FourierTransforms<stencil>;

}

}

}

}
