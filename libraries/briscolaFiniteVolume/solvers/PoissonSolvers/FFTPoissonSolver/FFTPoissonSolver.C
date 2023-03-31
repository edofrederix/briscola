#include "FFTPoissonSolver.H"
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
    cellSizes_(fvMsh.msh().cast<rectilinearMesh>().cellSizes()),
    decomp_(fvMsh),
    initData_(decomp_.Ni()[Pstream::myProcNo()]),
    xPencil_(decomp_.Nx()[Pstream::myProcNo()]),
    yPencil_(decomp_.Ny()[Pstream::myProcNo()]),
    zPencil_(decomp_.Nz()[Pstream::myProcNo()]),
    fft_(fvMsh, decomp_, xPencil_, yPencil_),
    tds_
    (
        fvMsh,
        decomp_.Nz()[Pstream::myProcNo()],
        decomp_.Sz()[Pstream::myProcNo()],
        fft_.BC()
    )
{
    checkMesh();
}

FFTPoissonSolver::FFTPoissonSolver
(
    const fvMesh& fvMsh
)
:
    PoissonSolver<stencil,scalar,colocated>(dictionary(),fvMsh),
    cellSizes_(fvMsh.msh().cast<rectilinearMesh>().cellSizes()),
    decomp_(fvMsh),
    initData_(decomp_.Ni()[Pstream::myProcNo()]),
    xPencil_(decomp_.Nx()[Pstream::myProcNo()]),
    yPencil_(decomp_.Ny()[Pstream::myProcNo()]),
    zPencil_(decomp_.Nz()[Pstream::myProcNo()]),
    fft_(fvMsh, decomp_, xPencil_, yPencil_),
    tds_
    (
        fvMsh,
        decomp_.Nz()[Pstream::myProcNo()],
        decomp_.Sz()[Pstream::myProcNo()],
        fft_.BC()
    )
{
    checkMesh();
}

void FFTPoissonSolver::solve
(
    colocatedScalarField& x,
    const colocatedScalarField* bPtr,
    const colocatedFaceScalarField* lambdaPtr,
    const bool ddt
)
{
    if (lambdaPtr != nullptr)
    {
        FatalErrorInFunction
            << "FFT Poisson solver does not support use of "
            << "lambda argument in solve function." << endl
            << abort(FatalError);
    }

    if (bPtr != nullptr)
    {
        // Initial decomposition
        labelVector I(decomp_.I());

        // Pencil decompositions
        labelVector X(decomp_.X());
        labelVector Y(decomp_.Y());
        labelVector Z(decomp_.Z());

        // Copy the data from bPtr to a scalarBlock
        // minus sign since bPtr = - RHS of the Poisson equation
        forAllCells((*bPtr)[0][0], i, j, k)
        {
            initData_(i,j,k) = - (*bPtr)[0][0](i,j,k);
        }

        // Transpose the RHS to x-pencils
        decomp_.transpose(initData_, xPencil_, I, X);

        // FFT of x-pencil data in x-direction
        fft_.fwdFFTx();

        // Transpose x-pencils to y-pencils
        decomp_.transpose(xPencil_, yPencil_, X, Y);

        // FFT of y-pencil data in y-direction
        fft_.fwdFFTy();

        // Transpose y-pencil data to z-pencils
        decomp_.transpose(yPencil_, zPencil_, Y, Z);

        // Solve tridiagonal systems in z-direction
        tds_.solve(zPencil_);

        // Transpose z-pencils to y-pencils
        decomp_.transpose(zPencil_, yPencil_, Z, Y);

        // Backward FFT of y-pencil data in y-direction
        fft_.bwdFFTy();

        // Transpose y-pencils to x-pencils
        decomp_.transpose(yPencil_, xPencil_, Y, X);

        // Backward FFT of x-pencil data in x-direction
        fft_.bwdFFTx();

        // Transpose x-pencils to initial decomposition
        decomp_.transpose(xPencil_, initData_, X, I);
    }
    else
    {
        initData_ *= 0;
    }

    // Copy scalarBlock values to pressure meshField
    forAllCells(x[0][0], i, j, k)
    {
        x[0][0](i,j,k) = initData_(i,j,k);
    }

    // Set ghost cells values
    x.correctBoundaryConditions();

    Info << "FFT: Solving for colocated p, residual = 0, nIter = 1" << endl;
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam
