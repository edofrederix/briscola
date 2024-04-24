#include "immersedBoundaryCondition.H"
#include "immersedBoundary.H"
#include "meshField.H"
#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"
#include "linearSystem.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

using Foam::max;
using Foam::min;
using Foam::sqr;
using Foam::mag;

// Constructor

template<class Type, class MeshType>
immersedBoundaryCondition<Type,MeshType>::immersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib,
    bool jac
)
:
    fvMsh_(mshField.fvMsh()),
    IB_(ib),
    dict_
    (
        mshField.found("boundaryConditions")
     && mshField.subDict("boundaryConditions").found(IB_.name())
      ? mshField.subDict("boundaryConditions").subDict(IB_.name())
      : dictionary::null
    ),
    JacobiGhostMethod_(jac),
    omega_
    (
        dict_.lookupOrDefault<scalar>("omega", 0.8)
    )
{}

template<class Type, class MeshType>
immersedBoundaryCondition<Type,MeshType>::~immersedBoundaryCondition()
{}

template<class Type, class MeshType>
autoPtr<immersedBoundaryCondition<Type,MeshType>>
immersedBoundaryCondition<Type,MeshType>::New
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
{
    dictionary IBCDict
    (
        mshField.subDict("boundaryConditions").subDict(ib.name())
    );

    const word IBCType(IBCDict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(IBCType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown immersed boundary condition "
            << IBCType << nl << nl
            << "Valid immersed boundary condition types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<immersedBoundaryCondition<Type,MeshType>>(cstrIter()(mshField, ib));
}

template<class Type, class MeshType>
void immersedBoundaryCondition<Type,MeshType>::correctLinearSystem
(
    linearSystem<stencil,Type,MeshType>&
)
{
    // Do nothing by default
}

template<class Type, class MeshType>
void immersedBoundaryCondition<Type,MeshType>::correctJacobiPoints
(
    meshLevel<Type,MeshType>&
) const
{
    // Do nothing by default
}

}

}

}
