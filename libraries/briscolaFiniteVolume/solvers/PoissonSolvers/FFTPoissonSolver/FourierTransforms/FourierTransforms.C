#include "FourierTransforms.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor
FourierTransforms::FourierTransforms
(
    const fvMesh& fvMsh,
    decomposer& d,
    scalarBlock& xPencil,
    scalarBlock& yPencil,
    scalarBlock& zPencil
)
:
    fvMsh_(fvMsh),
    N_(fvMsh.msh().cast<rectilinearMesh>().N()),
    decomp_(d),
    xPencil_(xPencil),
    yPencil_(yPencil),
    zPencil_(zPencil)
{
    FFTBoundaryConditions();
    FFTWplans();
}

// Destructor
FourierTransforms::~FourierTransforms()
{}

void FourierTransforms::FFTBoundaryConditions()
{

    colocatedScalarField p
    (
        "p",
        fvMsh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true,
        true
    );

    BC_ = labelVector(-1,-1,-1);

    switch (globalBoundaryConditionBaseType(p, faceOffsets[0]))
    {
        case 1:
            BC_.x() = 0;
            break;

        case 3:
            BC_.x() = 5; // C-C
            break;

        case 4:
            if (globalBoundaryConditionBaseType(p, faceOffsets[1]) == 4)
            {
                BC_.x() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[1]) == 5)
            {
                BC_.x() = 3; // D-N
            }
            break;

        case 5:
            if (globalBoundaryConditionBaseType(p, faceOffsets[1]) == 4)
            {
                BC_.x() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[1]) == 5)
            {
                BC_.x() = 2; // N-N
            }
            break;

        default:
            FatalError
                << "Incorrect pressure boundary condition." << endl
                << abort(FatalError);
            break;
    }

    switch (globalBoundaryConditionBaseType(p, faceOffsets[2]))
    {
        case 1:
            BC_.y() = 0;
            break;

        case 3:
            BC_.y() = 5; // C-C
            break;

        case 4:
            if (globalBoundaryConditionBaseType(p, faceOffsets[3]) == 4)
            {
                BC_.y() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[3]) == 5)
            {
                BC_.y() = 3; // D-N
            }
            break;

        case 5:
            if (globalBoundaryConditionBaseType(p, faceOffsets[3]) == 4)
            {
                BC_.y() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[3]) == 5)
            {
                BC_.y() = 2; // N-N
            }
            break;

        default:
            FatalError
                << "Incorrect pressure boundary condition." << endl
                << abort(FatalError);
            break;
    }

    switch (globalBoundaryConditionBaseType(p, faceOffsets[4]))
    {
        case 1:
            BC_.z() = 0;
            break;

        case 3:
            BC_.z() = 5; // C-C
            break;

        case 4:
            if (globalBoundaryConditionBaseType(p, faceOffsets[5]) == 4)
            {
                BC_.z() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[5]) == 5)
            {
                BC_.z() = 3; // D-N
            }
            break;

        case 5:
            if (globalBoundaryConditionBaseType(p, faceOffsets[5]) == 4)
            {
                BC_.z() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[5]) == 5)
            {
                BC_.z() = 2; // N-N
            }
            break;

        default:
            FatalError
                << "Incorrect pressure boundary condition." << endl
                << abort(FatalError);
            break;
    }
}

void FourierTransforms::FFTWplans()
{
    fftw_r2r_kind transforms[] =
    {
        // forward transform types
        FFTW_RODFT10, FFTW_RODFT10, FFTW_REDFT10, FFTW_RODFT11, FFTW_REDFT11, FFTW_R2HC,
        // backward transform types
        FFTW_RODFT01, FFTW_RODFT01, FFTW_REDFT01, FFTW_RODFT11, FFTW_REDFT11, FFTW_HC2R
    };

    // Transform type based on boundary conditions

    fftw_r2r_kind kind_fwd_x[] = {transforms[BC_.x()]};
    fftw_r2r_kind kind_bwd_x[] = {transforms[BC_.x()+6]};

    fftw_r2r_kind kind_fwd_y[] = {transforms[BC_.y()]};
    fftw_r2r_kind kind_bwd_y[] = {transforms[BC_.y()+6]};

    fftw_r2r_kind kind_fwd_z[] = {transforms[BC_.z()]};
    fftw_r2r_kind kind_bwd_z[] = {transforms[BC_.z()+6]};

    // Pencil dimensions

    labelVector Nx = decomp_.Nx()[Pstream::myProcNo()];
    labelVector Ny = decomp_.Ny()[Pstream::myProcNo()];
    labelVector Nz = decomp_.Nz()[Pstream::myProcNo()];

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
        reinterpret_cast<double*>(xPencil_.begin()),
        Na_x,
        stride,
        sep,
        reinterpret_cast<double*>(xPencil_.begin()),
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
        reinterpret_cast<double*>(xPencil_.begin()),
        Na_x,
        stride,
        sep,
        reinterpret_cast<double*>(xPencil_.begin()),
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
        reinterpret_cast<double*>(yPencil_.begin()),
        Na_y,
        stride,
        sep,
        reinterpret_cast<double*>(yPencil_.begin()),
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
        reinterpret_cast<double*>(yPencil_.begin()),
        Na_y,
        stride,
        sep,
        reinterpret_cast<double*>(yPencil_.begin()),
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
        reinterpret_cast<double*>(zPencil_.begin()),
        Na_z,
        stride,
        sep,
        reinterpret_cast<double*>(zPencil_.begin()),
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
        reinterpret_cast<double*>(zPencil_.begin()),
        Na_z,
        stride,
        sep,
        reinterpret_cast<double*>(zPencil_.begin()),
        Na_z,
        stride,
        sep,
        kind_bwd_z,
        FFTW_MEASURE
    );
}


void FourierTransforms::fwdFFTx()
{
    if (BC_.x()!= 0)
    {
        normalization_ = 1.0;
        fftw_execute(fwdPlanX_);
    }
}

void FourierTransforms::bwdFFTx()
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

void FourierTransforms::fwdFFTy()
{
    if (BC_.y() != 0)
    {
        normalization_ = 1.0;
        fftw_execute(fwdPlanY_);
    }
}

void FourierTransforms::bwdFFTy()
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

void FourierTransforms::fwdFFTz()
{
    if (BC_.z() != 0)
    {
        normalization_ = 1.0;
        fftw_execute(fwdPlanZ_);
    }
}

void FourierTransforms::bwdFFTz()
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

void FourierTransforms::normalize(scalarBlock& transformedData)
{
    transformedData /= normalization_;
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam