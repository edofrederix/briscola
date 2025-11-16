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
    {
        if (!sys.diagonal())
            FatalErrorInFunction
                << "System is not diagonal" << endl
                << abort(FatalError);
    }

    diagonalSmoother<SType,Type,MeshType>::Smooth
    (
        sys,
        0,
        1,
        labelList(MeshType::numberOfDirections, 0)
    );

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
