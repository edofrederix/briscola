#include "immersedBoundaryCondition.H"
#include "immersedBoundary.H"
#include "meshField.H"
#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"
#include "linearSystem.H"
#include "emptyImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
List<Type> immersedBoundaryCondition<Type,MeshType>::read(const word entry)
const
{
    typedef typename outerProduct<Type, typename MeshType::dimType>::type
        physicalType;

    const physicalType value = this->dict_.lookup<physicalType>(entry);

    List<Type> values(MeshType::numberOfDirections);

    tensor base(eye);

    if (this->fvMsh_.msh().template castable<rectilinearMesh>())
        base = this->fvMsh_.msh().template cast<rectilinearMesh>().base();

    for (int d = 0; d < MeshType::numberOfDirections; d++)
        values[d] = MeshType::project(value, d, base);

    return values;
}

// Constructor

template<class Type, class MeshType>
immersedBoundaryCondition<Type,MeshType>::immersedBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const immersedBoundary<MeshType>& ib,
    const meshField<label,MeshType>* maskPtr
)
:
    fvMsh_(field.fvMsh()),
    field_(const_cast<meshField<Type,MeshType>&>(field)),
    ib_(ib),
    forcingMaskPtr_(maskPtr),
    dict_
    (
        field.found("boundaryConditions")
     && field.subDict("boundaryConditions").found(ib_.name())
      ? field.subDict("boundaryConditions").subDict(ib_.name())
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
    const meshField<Type,MeshType>& field,
    const immersedBoundary<MeshType>& ib
)
{
    // For fields with an existing boundary conditions dictionary, the immersed
    // boundary condition must be specified. For other fields an empty immersed
    // boundary condition is returned.

    if (field.found("boundaryConditions"))
    {
        dictionary ibmDict
        (
            field.subDict("boundaryConditions").subDict(ib.name())
        );

        const word ibcType(ibmDict.lookup("type"));

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(ibcType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown immersed boundary condition "
                << ibcType << nl << nl
                << "Valid immersed boundary condition types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<immersedBoundaryCondition<Type,MeshType>>
        (
            cstrIter()(field, ib)
        );
    }
    else
    {
        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find
            (
                emptyImmersedBoundaryCondition<Type,MeshType>::typeName
            );

        return autoPtr<immersedBoundaryCondition<Type,MeshType>>
        (
            cstrIter()(field, ib)
        );
    }
}

template<class Type, class MeshType>
void immersedBoundaryCondition<Type,MeshType>::evaluate(const label l)
{
    forAll(field_[l], d)
        this->evaluate(l,d);
}

}

}

}
