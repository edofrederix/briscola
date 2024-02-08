#include "FFTPoissonSolver.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

typedef FFTPoissonSolver<stencil> FFTPoissonSolverStencil;
typedef FFTPoissonSolver<symmStencil> FFTPoissonSolverSymmStencil;

defineTemplateTypeNameAndDebugWithName(FFTPoissonSolverStencil,"FFT",0)
defineTemplateTypeNameAndDebugWithName(FFTPoissonSolverSymmStencil,"FFT",0)

PoissonSolver<stencil,scalar,colocated>::
adddictionaryConstructorToTable<FFTPoissonSolverStencil>
    addFFTPoissonSolverStencilConstructorToTable_;

PoissonSolver<symmStencil,scalar,colocated>::
adddictionaryConstructorToTable<FFTPoissonSolverSymmStencil>
    addFFTPoissonSolverSymmStencilConstructorToTable_;

template<class SType>
void FFTPoissonSolver<SType>::checkMesh(const fvMesh& fvMsh)
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

template<class SType>
FFTPoissonSolver<SType>::FFTPoissonSolver
(
    const word PoissonSolverName,
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    PoissonSolver<SType,scalar,colocated>(dict,fvMsh),
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

template<class SType>
FFTPoissonSolver<SType>::FFTPoissonSolver
(
    const fvMesh& fvMsh
)
:
    PoissonSolver<SType,scalar,colocated>(dictionary(),fvMsh),
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

template<class SType>
void FFTPoissonSolver<SType>::solve
(
    colocatedScalarField& x,
    const colocatedScalarField* bPtr,
    const colocatedLowerFaceScalarField* lambdaPtr,
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
        fft_.reset(new FFT::FourierTransforms<SType>(*this, x));
        tds_.reset(new FFT::tridiagonalSolver<SType>(*this, fft_->BC()));
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
    forAllCells(x, i, j, k)
    {
        if (bPtr != nullptr)
        {
            initData_(i,j,k) = - (*bPtr)(i,j,k);
        }

        if (ddt)
        {
            initData_(i,j,k) -= x.oldTime()(i,j,k)/deltaT;
        }
    }

    const labelVector N(fvMsh_.msh().cast<rectilinearMesh>().N());

    const PtrList<boundaryCondition<scalar,colocated>>& bcs =
        x.boundaryConditions();

    // Correct boundaries to evaluate boundary values
    x.correctBoundaryConditions();

    // Correct RHS for inhomogeneous BC's
    forAll(bcs, bci)
    {
        if (bcs[bci].offsetDegree() == 1)
        {
            const labelVector bo = bcs[bci].offset();

            const PtrList<PartialList<scalar>>& cellSizes =
                fvMsh_.msh().cast<rectilinearMesh>().globalCellSizes();

            const labelVector S(fvMsh_.template S<colocated>(bo));
            const labelVector E(fvMsh_.template E<colocated>(bo));

            labelVector ijk;

            const label dir(abs(bo.y())*1+abs(bo.z())*2);
            const label cellIndex((bo[dir] < 0) ? 0 : (N[dir]-1));
            const scalar sqrCellSize
            (
                Foam::sqr(cellSizes[dir][cellIndex])
            );

            if (bcs[bci].baseType() == DIRICHLETBC)
            {
                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    initData_(ijk)
                        -= (x(ijk+bo) + x(ijk))/sqrCellSize;
                }
            }

            if (bcs[bci].baseType() == NEUMANNBC)
            {
                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    initData_(ijk)
                        -= (x(ijk+bo) - x(ijk))/sqrCellSize;
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
    forAllCells(x, i, j, k)
    {
        x(i,j,k) = initData_(i,j,k);
    }

    // Set ghost cells values
    x.correctBoundaryConditions();

    Info<< "FFT: Solving for colocated " << x.name()
        << ", residual = 0, nIter = 1" << endl;
}

// Instantiate

template class FFTPoissonSolver<stencil>;
template class FFTPoissonSolver<symmStencil>;

}

}

}