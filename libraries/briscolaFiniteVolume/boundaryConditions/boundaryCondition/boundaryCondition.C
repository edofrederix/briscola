#include "boundaryCondition.H"

#include "colocated.H"
#include "staggered.H"

#include "stencil.H"
#include "symmStencil.H"
#include "diagStencil.H"

#include "meshField.H"

namespace Foam
{

template<>
const char* NamedEnum<briscola::fv::boundaryConditionBaseType,7>::names[] =
{
    "dummy",
    "empty",
    "parallel",
    "periodic",
    "Dirichlet",
    "Neumann",
    "Robin"
};

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
boundaryCondition<Type,MeshType>::boundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
:
    refCount(),
    fvMsh_(mshField.fvMsh()),
    mshField_(const_cast<meshField<Type,MeshType>&>(mshField)),
    patch_(patch),
    dict_
    (
        mshField_.found("boundaryConditions")
     && mshField_.subDict("boundaryConditions").found(patch_.name())
      ? mshField_.subDict("boundaryConditions").subDict(patch_.name())
      : dictionary::null
    )
{}


template<class Type, class MeshType>
boundaryCondition<Type,MeshType>::boundaryCondition
(
    const boundaryCondition<Type,MeshType>& bc
)
:
    refCount(),
    fvMsh_(bc.fvMsh_),
    mshField_(bc.mshField_),
    patch_(bc.patch_),
    dict_(bc.dict_)
{}

template<class Type, class MeshType>
boundaryCondition<Type,MeshType>::boundaryCondition
(
    const meshField<Type,MeshType>& field,
    const boundaryCondition<Type,MeshType>& bc
)
:
    refCount(),
    fvMsh_(bc.fvMsh_),
    mshField_(const_cast<meshField<Type,MeshType>&>(field)),
    patch_(bc.patch_),
    dict_(bc.dict_)
{}

template<class Type, class MeshType>
boundaryCondition<Type,MeshType>::~boundaryCondition()
{}

template<class Type, class MeshType>
autoPtr<boundaryCondition<Type,MeshType>> boundaryCondition<Type,MeshType>::New
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch,
    const word boundaryConditionType
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(boundaryConditionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown boundary condition type "
            << boundaryConditionType << nl << nl
            << "Valid boundary condition types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<boundaryCondition<Type,MeshType>>
    (
        cstrIter()(mshField, patch)
    );
}

template<class Type, class MeshType>
autoPtr<boundaryCondition<Type,MeshType>> boundaryCondition<Type,MeshType>::NewBoundary
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
{
    const word boundaryConditionType
    (
        mshField.found("boundaryConditions")
      ? word
        (
            mshField.subDict("boundaryConditions")
           .subDict(patch.name()).lookup("type")
        )
      : word("dummy")
    );

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(boundaryConditionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown boundary condition type "
            << boundaryConditionType << nl << nl
            << "Valid boundary condition types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<boundaryCondition<Type,MeshType>>
    (
        cstrIter()(mshField, patch)
    );
}

template<class Type, class MeshType>
autoPtr<boundaryCondition<Type,MeshType>>
boundaryCondition<Type,MeshType>::NewPeriodic
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find("periodic");

    return autoPtr<boundaryCondition<Type,MeshType>>
    (
        cstrIter()(mshField, patch)
    );
}

template<class Type, class MeshType>
autoPtr<boundaryCondition<Type,MeshType>>
boundaryCondition<Type,MeshType>::NewParallel
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find("parallel");

    return autoPtr<boundaryCondition<Type,MeshType>>
    (
        cstrIter()(mshField, patch)
    );
}

template<class Type, class MeshType>
const meshField<hexVector,MeshType>&
boundaryCondition<Type,MeshType>::faceCenters() const
{
    return fvMsh_.template metrics<MeshType>().faceCenters();
}

template<class Type, class MeshType>
const meshField<hexVector,MeshType>&
boundaryCondition<Type,MeshType>::faceNormals() const
{
    return fvMsh_.template metrics<MeshType>().faceNormals();
}

template<class Type, class MeshType>
const meshField<hexScalar,MeshType>&
boundaryCondition<Type,MeshType>::faceAreas() const
{
    return fvMsh_.template metrics<MeshType>().faceAreas();
}

template<class Type, class MeshType>
const meshField<scalar,MeshType>&
boundaryCondition<Type,MeshType>::cellVolumes() const
{
    return fvMsh_.template metrics<MeshType>().cellVolumes();
}

template<class Type, class MeshType>
const meshField<vector,MeshType>&
boundaryCondition<Type,MeshType>::cellCenters() const
{
    return fvMsh_.template metrics<MeshType>().cellCenters();
}

template<class Type, class MeshType>
const meshField<hexScalar,MeshType>&
boundaryCondition<Type,MeshType>::faceDeltas() const
{
    return fvMsh_.template metrics<MeshType>().faceDeltas();
}

// All meshField types must have the boundaryCondition class compiled. Most
// types don't actually have boundary conditions implemented.

makeBoundaryCondition(label,colocated)
makeBoundaryCondition(label,staggered)

makeBoundaryCondition(scalar,colocated)
makeBoundaryCondition(scalar,staggered)

makeBoundaryCondition(hexScalar,colocated)
makeBoundaryCondition(hexScalar,staggered)

makeBoundaryCondition(vector,colocated)
makeBoundaryCondition(vector,staggered)

makeBoundaryCondition(hexVector,colocated)
makeBoundaryCondition(hexVector,staggered)

makeBoundaryCondition(tensor,colocated)
makeBoundaryCondition(tensor,staggered)

makeBoundaryCondition(sphericalTensor,colocated)
makeBoundaryCondition(sphericalTensor,staggered)

makeBoundaryCondition(symmTensor,colocated)
makeBoundaryCondition(symmTensor,staggered)

makeBoundaryCondition(diagTensor,colocated)
makeBoundaryCondition(diagTensor,staggered)

makeBoundaryCondition(symmStencil,colocated)
makeBoundaryCondition(symmStencil,staggered)

makeBoundaryCondition(diagStencil,colocated)
makeBoundaryCondition(diagStencil,staggered)

makeBoundaryCondition(stencil,colocated)
makeBoundaryCondition(stencil,staggered)

}

}

}
