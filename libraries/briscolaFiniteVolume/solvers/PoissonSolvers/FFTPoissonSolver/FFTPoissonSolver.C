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
    PoissonSolver<stencil,scalar,colocated>(dict,fvMsh)
{
    checkMesh();
}

FFTPoissonSolver::FFTPoissonSolver
(
    const fvMesh& fvMsh
)
:
    PoissonSolver<stencil,scalar,colocated>(dictionary(),fvMsh)
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
    // FFT solver magic
}

}

}

}
