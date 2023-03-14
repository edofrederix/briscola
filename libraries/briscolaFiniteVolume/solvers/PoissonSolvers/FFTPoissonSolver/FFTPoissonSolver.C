#include "FFTPoissonSolver.H"
#include "imSchemes.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

void FFTPoissonSolver::checkMesh() const
{
    // Cast fails if the mesh is not rectilinear

    const rectilinearMesh& mesh = this->fvMsh_.msh().cast<rectilinearMesh>();

    if (cmptSum(mesh.uniform()) < 2)
    {
        FatalErrorInFunction
            << "At least two mesh directions must be uniform "
            << "for the " << this->type() << " solver." << endl
            << abort(FatalError);
    }
}

FFTPoissonSolver::FFTPoissonSolver
(
    const word PoissonSolverName,
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    PoissonSolver<stencil,scalar,colocated>(dict,fvMsh),
    cellSizes_ (fvMsh.msh().cast<rectilinearMesh>().cellSizes())
{
    checkMesh();

    decomp_ = new decomposer(fvMsh);

    xPencil_ = scalarBlock(decomp_->Nx()[Pstream::myProcNo()]);
    yPencil_ = scalarBlock(decomp_->Ny()[Pstream::myProcNo()]);

    fft_ = new FourierTransforms(fvMsh, *decomp_, xPencil_, yPencil_);

    tds_ = new tridiagonalSolver
    (
        fvMsh,
        decomp_->Nz()[Pstream::myProcNo()],
        decomp_->sz()[Pstream::myProcNo()],
        fft_->BC()
    );
}

FFTPoissonSolver::FFTPoissonSolver
(
    const fvMesh& fvMsh
)
:
    PoissonSolver<stencil,scalar,colocated>(dictionary(),fvMsh)
{
    checkMesh();


    decomp_ = new decomposer(fvMsh);

    xPencil_ = scalarBlock(decomp_->Nx()[Pstream::myProcNo()]);
    yPencil_ = scalarBlock(decomp_->Ny()[Pstream::myProcNo()]);

    fft_ = new FourierTransforms(fvMsh, *decomp_, xPencil_, yPencil_);

    tds_ = new tridiagonalSolver
    (
        fvMsh,
        decomp_->Nz()[Pstream::myProcNo()],
        decomp_->sz()[Pstream::myProcNo()],
        fft_->BC()
    );
}

void FFTPoissonSolver::solve
(
    colocatedScalarField& x,
    const colocatedScalarField* bPtr,
    const colocatedFaceScalarField* lambdaPtr,
    const bool ddt
)
{
    int rank = Pstream::myProcNo();

    // Mesh dimensions
    labelVector N(fvMsh_.msh().cast<rectilinearMesh>().N());

    // Global boundary conditions
    labelVector BC(fft_->BC());

    // Initial decomposition
    labelVector I(decomp_->I());
    List<labelVector> Ni = decomp_->Ni();
    List<labelVector> si = decomp_->si();

    // Pencil decompositions
    labelVector X(decomp_->X());
    List<labelVector> Nx = decomp_->Nx();
    List<labelVector> sx = decomp_->sx();

    labelVector Y(decomp_->Y());
    List<labelVector> Ny = decomp_->Ny();
    List<labelVector> sy = decomp_->sy();

    labelVector Z(decomp_->Z());
    List<labelVector> Nz = decomp_->Nz();
    List<labelVector> sz = decomp_->sz();
    scalarBlock zPencil(Nz[rank]);

    // Copy the data from bPtr to a scalarBlock
    // minus sign since bPtr = - RHS of the Poisson equation
    scalarBlock initData(Ni[rank]);

    for (int i = 0; i < Ni[rank].x(); i++)
    {
        for (int j = 0; j < Ni[rank].y(); j++)
        {
            for (int k = 0; k < Ni[rank].z(); k++)
            {
                initData(i,j,k) = - bPtr[0][0][0](i,j,k);
            }
        }
    }

    // Transpose the RHS to x-pencils
    decomp_->transpose(initData, xPencil_, I, X);

    // FFT of x-pencil data in x-direction
    fft_->fwdFFTx();

    // Transpose x-pencils to y-pencils
    decomp_->transpose(xPencil_, yPencil_, X, Y);

    // FFT of y-pencil data in y-direction
    decomp_->yTransFwd(yPencil_);

    fft_->fwdFFTy();

    decomp_->yTransBwd(yPencil_);

    // Transpose y-pencil data to z-pencils
    decomp_->transpose(yPencil_, zPencil, Y, Z);

    // Solve tridiagonal systems in z-direction
    tds_->solve(zPencil);

    // Transpose z-pencils to y-pencils
    decomp_->transpose(zPencil, yPencil_, Z, Y);

    // Backward FFT of y-pencil data in y-direction
    decomp_->yTransFwd(yPencil_);

    fft_->bwdFFTy();

    decomp_->yTransBwd(yPencil_);

    // Transpose y-pencils to x-pencils
    decomp_->transpose(yPencil_, xPencil_, Y, X);

    // Backward FFT of x-pencil data in x-direction
    fft_->bwdFFTx();

    // Transpose x-pencils to initial decomposition
    decomp_->transpose(xPencil_, initData, X, I);

    // Copy scalarBlock values to pressure meshField
    for (int i = 0; i < Ni[rank].x(); i++)
    {
        for (int j = 0; j < Ni[rank].y(); j++)
        {
            for (int k = 0; k < Ni[rank].z(); k++)
            {
                x[0][0](i,j,k) = initData(i,j,k);
            }
        }
    }

    // Set ghost cells values
    x.correctBoundaryConditions();

    Info << "FFT: Solving for colocated p, residual = 0, nIter = 1" << endl;
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam
