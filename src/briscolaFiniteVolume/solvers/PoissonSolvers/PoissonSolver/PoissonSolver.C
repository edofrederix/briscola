#include "PoissonSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
void PoissonSolver<SType,Type,MeshType>::initFlux()
{
    if (fluxPtr_.empty())
        fluxPtr_.reset
        (
            new faceField<Type,MeshType>
            (
                "PoissonFlux",
                fvMsh_
            )
        );
}

template<class SType, class Type, class MeshType>
PoissonSolver<SType,Type,MeshType>::PoissonSolver
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    dict_(dict),
    fvMsh_(fvMsh),
    computeFlux_(false),
    rkSchemePtr_(nullptr)
{
    if (fvMsh.db().foundObject<RungeKuttaScheme>("rkScheme"))
        rkSchemePtr_ =
            &fvMsh.db().lookupObjectRef<RungeKuttaScheme>("rkScheme");
}

template<class SType, class Type, class MeshType>
PoissonSolver<SType,Type,MeshType>::~PoissonSolver()
{}

template<class SType, class Type, class MeshType>
autoPtr<PoissonSolver<SType,Type,MeshType>>
PoissonSolver<SType,Type,MeshType>::New
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
{
    const word PoissonSolverType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(PoissonSolverType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        cstrIter = dictionaryConstructorTablePtr_->find("default");
    }

    return autoPtr<PoissonSolver<SType,Type,MeshType>>
    (
        cstrIter()(dict, fvMsh)
    );
}

template<class SType, class Type, class MeshType>
autoPtr<PoissonSolver<SType,Type,MeshType>> PoissonSolver<SType,Type,MeshType>::New
(
    const word PoissonSolverName,
    const fvMesh& fvMsh
)
{
    const dictionary dict
    (
        fvMsh.solverDict().subDict(PoissonSolverName)
    );

    return New(dict, fvMsh);
}

}

}

}
