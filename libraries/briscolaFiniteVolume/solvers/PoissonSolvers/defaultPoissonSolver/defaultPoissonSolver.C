#include "defaultPoissonSolver.H"
#include "imSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
defaultPoissonSolver<SType,Type,MeshType>::defaultPoissonSolver
(
    const word PoissonSolverName,
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    PoissonSolver<SType,Type,MeshType>(dict,fvMsh),
    solverPtr_
    (
        word(dict.lookup("type")) == "default"
      ? solver<SType,Type,MeshType>::New
        (
            dict.subDict("solver"),
            fvMsh
        )
      : solver<SType,Type,MeshType>::New
        (
            PoissonSolverName,
            fvMsh
        )
    )
{}

template<class SType, class Type, class MeshType>
void defaultPoissonSolver<SType,Type,MeshType>::solve
(
    meshField<Type,MeshType>& x,
    const meshField<Type,MeshType>* bPtr,
    const meshField<faceScalar,MeshType>* lambdaPtr,
    const bool ddt
)
{
    linearSystem<SType,Type,MeshType> sys(x);

    if (lambdaPtr)
    {
        sys = im::laplacian(*lambdaPtr,x);
    }
    else
    {
        sys = im::laplacian(x);
    }

    if (bPtr)
    {
        sys += (*bPtr);
    }

    if (ddt)
    {
        sys -= im::ddt(x);
    }

    solverPtr_->solve(sys);
}

}

}

}
