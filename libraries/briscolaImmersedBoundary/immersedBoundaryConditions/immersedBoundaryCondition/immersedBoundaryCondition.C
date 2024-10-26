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

// Constructor

template<class Type, class MeshType>
immersedBoundaryCondition<Type,MeshType>::immersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib,
    const meshField<label,MeshType>* maskPtr
)
:
    fvMsh_(mshField.fvMsh()),
    mshField_(const_cast<meshField<Type,MeshType>&>(mshField)),
    ib_(ib),
    forcingMaskPtr_(maskPtr),
    dict_
    (
        mshField.found("boundaryConditions")
     && mshField.subDict("boundaryConditions").found(ib_.name())
      ? mshField.subDict("boundaryConditions").subDict(ib_.name())
      : dictionary::null
    ),
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

    return autoPtr<immersedBoundaryCondition<Type,MeshType>>
    (
        cstrIter()(mshField, ib)
    );
}

template<class Type, class MeshType>
void immersedBoundaryCondition<Type,MeshType>::evaluate(const label l)
{
    forAll(mshField_[l], d)
        this->evaluate(l,d);
}

}

}

}
