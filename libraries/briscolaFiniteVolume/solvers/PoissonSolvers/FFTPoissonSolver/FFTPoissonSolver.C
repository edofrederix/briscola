#include "FFTPoissonSolver.H"
#include "imSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

void FFTPoissonSolver::checkMesh() const
{
    if (this->fvMsh_.rectilinear() != unitXYZ)
    {
        FatalErrorInFunction
            << "Mesh must be rectilinear in three directions "
            << "for the " << this->type() << " solver." << endl
            << abort(FatalError);
    }

    if (cmptSum(this->fvMsh_.uniform()) < 2)
    {
        FatalErrorInFunction
            << "At least two mesh directions must be uniform "
            << "for the " << this->type() << " solver." << endl
            << abort(FatalError);
    }
}

void FFTPoissonSolver::prepare()
{
    cellSizes_.clear();
    cellSizes_.setSize(3);

    for (int dir = 0; dir < 3; dir++)
    {
        cellSizes_.set
        (
            dir,
            new scalarList(this->fvMsh_.rectilinearCellSizes(dir))
        );
    }
}

FFTPoissonSolver::FFTPoissonSolver
(
    const word PoissonSolverName,
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    PoissonSolver<stencil,scalar,colocated>(dict,fvMsh)
{
    checkMesh();
    prepare();
}

FFTPoissonSolver::FFTPoissonSolver
(
    const fvMesh& fvMsh
)
:
    PoissonSolver<stencil,scalar,colocated>(dictionary(),fvMsh)
{
    checkMesh();
    prepare();
}

void FFTPoissonSolver::solve
(
    colocatedScalarField& x,
    const colocatedScalarField* bPtr,
    const colocatedFaceScalarField* lambdaPtr,
    const bool ddt
)
{
    // FFT solver magic
}

}

}

}
