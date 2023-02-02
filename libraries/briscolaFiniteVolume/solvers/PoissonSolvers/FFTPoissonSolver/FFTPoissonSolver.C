#include "FFTPoissonSolver.H"
#include "imSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
FFTPoissonSolver<SType,Type,MeshType>::FFTPoissonSolver
(
    const word PoissonSolverName,
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    PoissonSolver<SType,Type,MeshType>(dict,fvMsh)
{}

template<class SType, class Type, class MeshType>
void FFTPoissonSolver<SType,Type,MeshType>::solve
(
    meshField<Type,MeshType>& x,
    const meshField<Type,MeshType>* bPtr,
    const meshField<faceScalar,MeshType>* lambdaPtr,
    const bool ddt
)
{
    // FFT solver magic
}

}

}

}
