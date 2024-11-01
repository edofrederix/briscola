#include "diagonal.H"

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
    sys.eliminateGhosts();
    sys.setForcingMask();

    sys.x().makeShallow();
    sys.b().makeShallow();

    List<bool> diag = sys.diagonal();

    forAll(diag, i)
        if (!diag[i])
            FatalErrorInFunction
                << "Direction " << i << " of " << sys.name()
                << " is not diagonal." << endl
                << abort(FatalError);

    for (int d = 0; d < MeshType::numberOfDirections; d++)
        solver<SType,Type,MeshType>::smoother::smoothDiag(sys, 0, d);

    sys.x().correctBoundaryConditions();

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
