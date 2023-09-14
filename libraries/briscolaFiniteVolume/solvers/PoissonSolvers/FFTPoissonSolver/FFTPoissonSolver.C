#include "FFTPoissonSolver.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(FFTPoissonSolver,0);

PoissonSolver<stencil,scalar,colocated>::
adddictionaryConstructorToTable<FFTPoissonSolver>
    addFFTPoissonSolverConstructorToTable_;

void FFTPoissonSolver::checkMesh(const fvMesh& fvMsh)
{
    const rectilinearMesh& mesh = fvMsh.msh().cast<rectilinearMesh>();

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
    fvMsh_(fvMsh),
    FFTPlan_(fvMsh),
    decomp_(fvMsh, FFTPlan_.decompType()),
    initData_(decomp_.Ni()[Pstream::myProcNo()]),
    xPencil_(decomp_.Nx()[Pstream::myProcNo()]),
    yPencil_(decomp_.Ny()[Pstream::myProcNo()]),
    zPencil_(decomp_.Nz()[Pstream::myProcNo()])
{
    checkMesh(fvMsh);
}

FFTPoissonSolver::FFTPoissonSolver
(
    const fvMesh& fvMsh
)
:
    PoissonSolver<stencil,scalar,colocated>(dictionary(),fvMsh),
    fvMsh_(fvMsh),
    FFTPlan_(fvMsh),
    decomp_(fvMsh, FFTPlan_.decompType()),
    initData_(decomp_.Ni()[Pstream::myProcNo()], Zero),
    xPencil_(decomp_.Nx()[Pstream::myProcNo()]),
    yPencil_(decomp_.Ny()[Pstream::myProcNo()]),
    zPencil_(decomp_.Nz()[Pstream::myProcNo()])
{
    checkMesh(fvMsh);
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

    if (fft_.empty() || &fft_->x() != &x)
    {
        fft_.reset(new FFT::FourierTransforms(*this, x));
        tds_.reset(new FFT::tridiagonalSolver(*this, fft_->BC()));
    }
    // Initial decomposition
    labelVector I(decomp_.I());

    // Pencil decompositions
    labelVector X(decomp_.X());
    labelVector Y(decomp_.Y());
    labelVector Z(decomp_.Z());

    const scalar deltaT = x.fvMsh().time().deltaTValue();

    // Copy the data from bPtr to a scalarBlock
    // minus sign since bPtr = - RHS of the Poisson equation
    forAllCells(x[0][0], i, j, k)
    {
        if (bPtr != nullptr)
        {
            initData_(i,j,k) = - (*bPtr)[0][0](i,j,k);
        }

        if (ddt)
        {
            initData_(i,j,k) -= x.oldTime()[0][0](i,j,k)/deltaT;
        }
    }

    const labelVector N(fvMsh_.msh().cast<rectilinearMesh>().N());

    const PtrList<boundaryCondition<scalar,colocated>>& bcs =
        x.boundaryConditions();

    forAll(bcs, bci)
    {
        const labelVector bo = bcs[bci].boundaryOffset();

        if (bcs[bci].baseType() == 4)
        {
            // Get Dirichlet BC value
            scalarList values(bcs[bci].dict().lookup("values"));
            scalar value(values[0]);

            // Correct inhomogeneous BC
            if (value != 0)
            {
                const PtrList<PartialList<scalar>>& cellSizes =
                    fvMsh_.msh().cast<rectilinearMesh>().globalCellSizes();

                scalar cellSize = -1;
                if (bo.x() == -1)
                {
                    cellSize = cellSizes[0][0];

                    for (int j = 0; j < initData_.m(); j++)
                    for (int k = 0; k < initData_.n(); k++)
                    {
                        initData_(0,j,k) -=
                            2.0*value/Foam::sqr(cellSize);
                    }
                }
                else if (bo.x() == 1)
                {
                    cellSize = cellSizes[0][N.x()-1];

                    for (int j = 0; j < initData_.m(); j++)
                    for (int k = 0; k < initData_.n(); k++)
                    {
                        initData_(initData_.l()-1,j,k) -=
                            2.0*value/Foam::sqr(cellSize);
                    }
                }
                else if (bo.y() == -1)
                {
                    cellSize = cellSizes[1][0];

                    for (int i = 0; i < initData_.l(); i++)
                    for (int k = 0; k < initData_.n(); k++)
                    {
                        initData_(i,0,k) -=
                            2.0*value/Foam::sqr(cellSize);
                    }
                }
                else if (bo.y() == 1)
                {
                    cellSize = cellSizes[1][N.y()-1];

                    for (int i = 0; i < initData_.l(); i++)
                    for (int k = 0; k < initData_.n(); k++)
                    {
                        initData_(i,initData_.m()-1,k) -=
                            2.0*value/Foam::sqr(cellSize);
                    }
                }
                else if (bo.z() == -1)
                {
                    cellSize = cellSizes[2][0];

                    for (int i = 0; i < initData_.l(); i++)
                    for (int j = 0; j < initData_.m(); j++)
                    {
                        initData_(i,j,0) -=
                            2.0*value/Foam::sqr(cellSize);
                    }
                }
                else if (bo.z() == 1)
                {
                    cellSize = cellSizes[2][N.z()-1];

                    for (int i = 0; i < initData_.l(); i++)
                    for (int j = 0; j < initData_.m(); j++)
                    {
                        initData_(i,j,initData_.n()-1) -=
                            2.0*value/Foam::sqr(cellSize);
                    }
                }
            }
        }

        if (bcs[bci].baseType() == 5)
        {
            // Get Neumann BC gradient
            scalarList gradients(bcs[bci].dict().lookup("gradients"));
            scalar gradient(gradients[0]);

            // Correct inhomogeneous BC
            if (gradient != 0)
            {
                const PtrList<PartialList<scalar>>& cellSizes =
                    fvMsh_.msh().cast<rectilinearMesh>().globalCellSizes();

                scalar cellSize = -1;
                if (bo.x() == -1)
                {
                    cellSize = cellSizes[0][0];

                    for (int j = 0; j < initData_.m(); j++)
                    for (int k = 0; k < initData_.n(); k++)
                    {
                        initData_(0,j,k) -=
                            gradient/cellSize;
                    }
                }
                else if (bo.x() == 1)
                {
                    cellSize = cellSizes[0][N.x()-1];

                    for (int j = 0; j < initData_.m(); j++)
                    for (int k = 0; k < initData_.n(); k++)
                    {
                        initData_(initData_.l()-1,j,k) -=
                            gradient/cellSize;
                    }
                }
                else if (bo.y() == -1)
                {
                    cellSize = cellSizes[1][0];

                    for (int i = 0; i < initData_.l(); i++)
                    for (int k = 0; k < initData_.n(); k++)
                    {
                        initData_(i,0,k) -=
                            gradient/cellSize;
                    }
                }
                else if (bo.y() == 1)
                {
                    cellSize = cellSizes[1][N.y()-1];

                    for (int i = 0; i < initData_.l(); i++)
                    for (int k = 0; k < initData_.n(); k++)
                    {
                        initData_(i,initData_.m()-1,k) -=
                            gradient/cellSize;
                    }
                }
                else if (bo.z() == -1)
                {
                    cellSize = cellSizes[2][0];

                    for (int i = 0; i < initData_.l(); i++)
                    for (int j = 0; j < initData_.m(); j++)
                    {
                        initData_(i,j,0) -=
                            gradient/cellSize;
                    }
                }
                else if (bo.z() == 1)
                {
                    cellSize = cellSizes[2][N.z()-1];

                    for (int i = 0; i < initData_.l(); i++)
                    for (int j = 0; j < initData_.m(); j++)
                    {
                        initData_(i,j,initData_.n()-1) -=
                            gradient/cellSize;
                    }
                }
            }
        }
    }


    // Transform in two directions and solve in the third direction
    switch (FFTPlan_.solveDir())
    {
        case 0:
            if (FFTPlan_.firstTransDir() == 2)
            {
                decomp_.transpose(initData_, zPencil_, I, Z, "z", "z");

                fft_->fwdFFTz();

                if (N.y() > 1)
                {
                    decomp_.transpose(zPencil_, yPencil_, Z, Y, "z", "y");

                    fft_->fwdFFTy();

                    decomp_.transpose(yPencil_, xPencil_, Y, X, "y", "x");
                }
                else
                {
                    decomp_.transpose(zPencil_, xPencil_, Z, X, "z", "x");
                }

                tds_->solve(xPencil_, ddt);

                if (N.y() > 1)
                {
                    decomp_.transpose(xPencil_, yPencil_, X, Y, "x", "y");

                    fft_->bwdFFTy();

                    decomp_.transpose(yPencil_, zPencil_, Y, Z, "y", "z");
                }
                else
                {
                    decomp_.transpose(xPencil_, zPencil_, X, Z, "x", "z");
                }

                fft_->bwdFFTz();

                decomp_.transpose(zPencil_, initData_, Z, I, "z", "z");
            }
            else
            {
                decomp_.transpose(initData_, yPencil_, I, Y, "z", "y");

                fft_->fwdFFTy();

                if (N.z() > 1)
                {
                    decomp_.transpose(yPencil_, zPencil_, Y, Z, "y", "z");

                    fft_->fwdFFTz();

                    decomp_.transpose(zPencil_, xPencil_, Z, X, "z", "x");
                }
                else
                {
                    decomp_.transpose(yPencil_, xPencil_, Y, X, "y", "x");
                }

                tds_->solve(xPencil_, ddt);

                if (N.z() > 1)
                {
                    decomp_.transpose(xPencil_, zPencil_, X, Z, "x", "z");

                    fft_->bwdFFTz();

                    decomp_.transpose(zPencil_, yPencil_, Z, Y, "z", "y");
                }
                else
                {
                    decomp_.transpose(xPencil_, yPencil_, X, Y, "x", "y");
                }

                fft_->bwdFFTy();

                decomp_.transpose(yPencil_, initData_, Y, I, "y", "z");
            }
            break;

        case 1:
            if (FFTPlan_.firstTransDir() == 2)
            {
                decomp_.transpose(initData_, zPencil_, I, Z, "z", "z");

                fft_->fwdFFTz();

                if (N.x() > 1)
                {
                    decomp_.transpose(zPencil_, xPencil_, Z, X, "z", "x");

                    fft_->fwdFFTx();

                    decomp_.transpose(xPencil_, yPencil_, X, Y, "x", "y");
                }
                else
                {
                    decomp_.transpose(zPencil_, yPencil_, Z, Y, "z", "y");
                }

                tds_->solve(yPencil_, ddt);

                if (N.x() > 1)
                {
                    decomp_.transpose(yPencil_, xPencil_, Y, X, "y", "x");

                    fft_->bwdFFTx();

                    decomp_.transpose(xPencil_, zPencil_, X, Z, "x", "z");
                }
                else
                {
                    decomp_.transpose(yPencil_, zPencil_, Y, Z, "y", "z");
                }

                fft_->bwdFFTz();

                decomp_.transpose(zPencil_, initData_, Z, I, "z", "z");
            }
            else
            {
                decomp_.transpose(initData_, xPencil_, I, X, "z", "x");

                fft_->fwdFFTx();

                if (N.z() > 1)
                {
                    decomp_.transpose(xPencil_, zPencil_, X, Z, "x", "z");

                    fft_->fwdFFTz();

                    decomp_.transpose(zPencil_, yPencil_, Z, Y, "z", "y");
                }
                else
                {
                    decomp_.transpose(xPencil_, yPencil_, X, Y, "x", "y");
                }

                tds_->solve(yPencil_, ddt);

                if (N.z() > 1)
                {
                    decomp_.transpose(yPencil_, zPencil_, Y, Z, "y", "z");

                    fft_->bwdFFTz();

                    decomp_.transpose(zPencil_, xPencil_, Z, X, "z", "x");
                }
                else
                {
                    decomp_.transpose(yPencil_, xPencil_, Y, X, "y", "x");
                }

                fft_->bwdFFTx();

                decomp_.transpose(xPencil_, initData_, X, I, "x", "z");
            }
            break;

        case 2:
            if (FFTPlan_.firstTransDir() == 1)
            {
                decomp_.transpose(initData_, yPencil_, I, Y, "z", "y");

                fft_->fwdFFTy();

                if (N.x() > 1)
                {
                    decomp_.transpose(yPencil_, xPencil_, Y, X, "y", "x");

                    fft_->fwdFFTx();

                    decomp_.transpose(xPencil_, zPencil_, X, Z, "x", "z");
                }
                else
                {
                    decomp_.transpose(yPencil_, zPencil_, Y, Z, "y", "z");
                }

                tds_->solve(zPencil_, ddt);

                if (N.x() > 1)
                {
                    decomp_.transpose(zPencil_, xPencil_, Z, X, "z", "x");

                    fft_->bwdFFTx();

                    decomp_.transpose(xPencil_, yPencil_, X, Y, "x", "y");
                }
                else
                {
                    decomp_.transpose(zPencil_, yPencil_, Z, Y, "z", "y");
                }

                fft_->bwdFFTy();

                decomp_.transpose(yPencil_, initData_, Y, I, "y", "z");
            }
            else
            {
                decomp_.transpose(initData_, xPencil_, I, X, "z", "x");

                fft_->fwdFFTx();

                if (N.y() > 1)
                {
                    decomp_.transpose(xPencil_, yPencil_, X, Y, "x", "y");

                    fft_->fwdFFTy();

                    decomp_.transpose(yPencil_, zPencil_, Y, Z, "y", "z");
                }
                else
                {
                    decomp_.transpose(xPencil_, zPencil_, X, Z, "x", "z");
                }

                tds_->solve(zPencil_, ddt);


                if (N.y() > 1)
                {
                    decomp_.transpose(zPencil_, yPencil_, Z, Y, "z", "y");

                    fft_->bwdFFTy();

                    decomp_.transpose(yPencil_, xPencil_, Y, X, "y", "x");
                }
                else
                {
                    decomp_.transpose(zPencil_, xPencil_, Z, X, "z", "x");
                }

                fft_->bwdFFTx();

                decomp_.transpose(xPencil_, initData_, X, I, "x", "z");
            }
            break;

        default:
            FatalError
                << "Invalid solve direction."
                << endl;
            FatalError.exit();
            break;
    }

    fft_->normalize(initData_);

    // Copy scalarBlock values to solution meshField
    forAllCells(x[0][0], i, j, k)
    {
        x[0][0](i,j,k) = initData_(i,j,k);
    }

    // Set ghost cells values
    x.correctBoundaryConditions();

    Info<< "FFT: Solving for colocated " << x.name()
        << ", residual = 0, nIter = 1" << endl;
}

}

}

}