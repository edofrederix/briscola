#include "diagonal.H"
#include "diagonalSmoother.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
diagonal<SType,Type,MeshType>::diagonal
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>(dict,fvMsh)
{}

template<class SType, class Type, class MeshType>
void diagonal<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys,
    const bool constMatrix
)
{
    if (SType::nCsComponents > 1)
        sys.eliminateGhosts();

    sys.setForcingMask();

    if (SType::nCsComponents > 1)
        if (!sys.diagonal())
            FatalErrorInFunction
                << "System is not diagonal" << endl
                << abort(FatalError);

    for (int d = 0; d < MeshType::numberOfDirections; d++)
        diagonalSmoother<SType,Type,MeshType>::Sweep(sys, 0, d);

    sys.x()[0].correctBoundaryConditions();

    this->printSolverStats
    (
        this->typeName,
        sys.x().name(),
        List<Type>(MeshType::numberOfDirections, pTraits<Type>::one),
        List<Type>(MeshType::numberOfDirections, Zero),
        0
    );
}

}

}

}
