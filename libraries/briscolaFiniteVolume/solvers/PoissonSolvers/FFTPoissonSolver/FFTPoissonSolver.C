#include "FFTPoissonSolver.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

void FFTPoissonSolver::checkMesh()
{
    if (cmptSum(meshUniform_) < 2)
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
    PoissonPlan_(fvMsh),
    decompType_(PoissonPlan_.decompType()),
    meshUniform_(PoissonPlan_.meshUniform()),
    solveDir_(PoissonPlan_.solveDir()),
    decomp_(fvMsh, decompType_),
    initData_(decomp_.Ni()[Pstream::myProcNo()]),
    xPencil_(decomp_.Nx()[Pstream::myProcNo()]),
    yPencil_(decomp_.Ny()[Pstream::myProcNo()]),
    zPencil_(decomp_.Nz()[Pstream::myProcNo()]),
    fft_(fvMsh, decomp_, xPencil_, yPencil_, zPencil_),
    tds_
    (
        fvMsh,
        solveDir_,
        decomp_.Nd(solveDir_)[Pstream::myProcNo()],
        decomp_.Sd(solveDir_)[Pstream::myProcNo()],
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
    PoissonPlan_(fvMsh),
    decompType_(PoissonPlan_.decompType()),
    meshUniform_(PoissonPlan_.meshUniform()),
    solveDir_(PoissonPlan_.solveDir()),
    decomp_(fvMsh, decompType_),
    initData_(decomp_.Ni()[Pstream::myProcNo()]),
    xPencil_(decomp_.Nx()[Pstream::myProcNo()]),
    yPencil_(decomp_.Ny()[Pstream::myProcNo()]),
    zPencil_(decomp_.Nz()[Pstream::myProcNo()]),
    fft_(fvMsh, decomp_, xPencil_, yPencil_, zPencil_),
    tds_
    (
        fvMsh,
        solveDir_,
        decomp_.Nd(solveDir_)[Pstream::myProcNo()],
        decomp_.Sd(solveDir_)[Pstream::myProcNo()],
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

        // Transform in two directions and solve in the third direction
        switch (solveDir_)
        {
            case 0:
                if (decompType_ == 4 || decompType_ == 6)
                {
                    decomp_.transpose(initData_, zPencil_, I, Z, "z", "z");

                    fft_.fwdFFTz();

                    decomp_.transpose(zPencil_, yPencil_, Z, Y, "z", "y");

                    fft_.fwdFFTy();

                    decomp_.transpose(yPencil_, xPencil_, Y, X, "y", "x");

                    tds_.solve(xPencil_);

                    decomp_.transpose(xPencil_, yPencil_, X, Y, "x", "y");

                    fft_.bwdFFTy();

                    decomp_.transpose(yPencil_, zPencil_, Y, Z, "y", "z");

                    fft_.bwdFFTz();

                    decomp_.transpose(zPencil_, initData_, Z, I, "z", "z");
                }
                else
                {
                    decomp_.transpose(initData_, yPencil_, I, Y, "z", "y");

                    fft_.fwdFFTy();

                    decomp_.transpose(yPencil_, zPencil_, Y, Z, "y", "z");

                    fft_.fwdFFTz();

                    decomp_.transpose(zPencil_, xPencil_, Z, X, "z", "x");

                    tds_.solve(xPencil_);

                    decomp_.transpose(xPencil_, zPencil_, X, Z, "x", "z");

                    fft_.bwdFFTz();

                    decomp_.transpose(zPencil_, yPencil_, Z, Y, "z", "y");

                    fft_.bwdFFTy();

                    decomp_.transpose(yPencil_, initData_, Y, I, "y", "z");
                }
                break;

            case 1:
                if (decompType_ == 4 || decompType_ == 7)
                {
                    decomp_.transpose(initData_, zPencil_, I, Z, "z", "z");

                    fft_.fwdFFTz();

                    decomp_.transpose(zPencil_, xPencil_, Z, X, "z", "x");

                    fft_.fwdFFTx();

                    decomp_.transpose(xPencil_, yPencil_, X, Y, "x", "y");

                    tds_.solve(yPencil_);

                    decomp_.transpose(yPencil_, xPencil_, Y, X, "y", "x");

                    fft_.bwdFFTx();

                    decomp_.transpose(xPencil_, zPencil_, X, Z, "x", "z");

                    fft_.bwdFFTz();

                    decomp_.transpose(zPencil_, initData_, Z, I, "z", "z");
                }
                else
                {
                    decomp_.transpose(initData_, xPencil_, I, X, "z", "x");

                    fft_.fwdFFTx();

                    decomp_.transpose(xPencil_, zPencil_, X, Z, "x", "z");

                    fft_.fwdFFTz();

                    decomp_.transpose(zPencil_, yPencil_, Z, Y, "z", "y");

                    tds_.solve(yPencil_);

                    decomp_.transpose(yPencil_, zPencil_, Y, Z, "y", "z");

                    fft_.bwdFFTz();

                    decomp_.transpose(zPencil_, xPencil_, Z, X, "z", "x");

                    fft_.bwdFFTx();

                    decomp_.transpose(xPencil_, initData_, X, I, "x", "z");
                }
                break;

            case 2:
                if (decompType_ == 3 || decompType_ == 7)
                {
                    decomp_.transpose(initData_, yPencil_, I, Y, "z", "y");

                    fft_.fwdFFTy();

                    decomp_.transpose(yPencil_, xPencil_, Y, X, "y", "x");

                    fft_.fwdFFTx();

                    decomp_.transpose(xPencil_, zPencil_, X, Z, "x", "z");

                    tds_.solve(zPencil_);

                    decomp_.transpose(zPencil_, xPencil_, Z, X, "z", "x");

                    fft_.bwdFFTx();

                    decomp_.transpose(xPencil_, yPencil_, X, Y, "x", "y");

                    fft_.bwdFFTy();

                    decomp_.transpose(yPencil_, initData_, Y, I, "y", "z");
                }
                else
                {
                    decomp_.transpose(initData_, xPencil_, I, X, "z", "x");

                    fft_.fwdFFTx();

                    decomp_.transpose(xPencil_, yPencil_, X, Y, "x", "y");

                    fft_.fwdFFTy();

                    decomp_.transpose(yPencil_, zPencil_, Y, Z, "y", "z");

                    tds_.solve(zPencil_);

                    decomp_.transpose(zPencil_, yPencil_, Z, Y, "z", "y");

                    fft_.bwdFFTy();

                    decomp_.transpose(yPencil_, xPencil_, Y, X, "y", "x");

                    fft_.bwdFFTx();

                    decomp_.transpose(xPencil_, initData_, X, I, "x", "z");
                }
        }
    }
    else
    {
        initData_ *= 0;
    }

    fft_.normalize(initData_);

    // Copy scalarBlock values to pressure meshField
    forAllCells(x[0][0], i, j, k)
    {
        x[0][0](i,j,k) = initData_(i,j,k);
    }

    // Set ghost cells values
    x.correctBoundaryConditions();

    Info << "FFT: Solving for colocated p, residual = 0, nIter = 1" << endl;
}

}

}

}
