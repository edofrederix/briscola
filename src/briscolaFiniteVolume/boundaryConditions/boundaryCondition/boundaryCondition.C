#include "boundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshLevel.H"
#include "linearSystem.H"
#include "patchBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
boundaryCondition<Type,MeshType>::boundaryCondition
(
    const meshLevel<Type,MeshType>& level,
    const boundary& b
)
:
    refCount(),
    fvMsh_(level.fvMsh()),
    level_(const_cast<meshLevel<Type,MeshType>&>(level)),
    l_(level_.levelNum()),
    b_(b),
    dict_(dictionary::null)
{
    meshField<Type,MeshType>& field = level_.field();

    const word name(IOobject::member(b_.name()));

    if
    (
        field.found("boundaryConditions")
     && field.subDict("boundaryConditions").found(name)
    )
    {
        dict_ = field.subDict("boundaryConditions").subDict(name);
    }
}

template<class Type, class MeshType>
boundaryCondition<Type,MeshType>::boundaryCondition
(
    const boundaryCondition<Type,MeshType>& bc
)
:
    refCount(),
    fvMsh_(bc.fvMsh_),
    level_(bc.level_),
    l_(bc.l_),
    b_(bc.b_),
    dict_(bc.dict_)
{}

template<class Type, class MeshType>
boundaryCondition<Type,MeshType>::boundaryCondition
(
    const boundaryCondition<Type,MeshType>& bc,
    const meshLevel<Type,MeshType>& level
)
:
    refCount(),
    fvMsh_(bc.fvMsh_),
    level_(const_cast<meshLevel<Type,MeshType>&>(level)),
    l_(level.levelNum()),
    b_(level.lvl().boundaries().find(bc.b_.name())),
    dict_(bc.dict_)
{}

template<class Type, class MeshType>
boundaryCondition<Type,MeshType>::~boundaryCondition()
{}

template<class Type, class MeshType>
autoPtr<boundaryCondition<Type,MeshType>> boundaryCondition<Type,MeshType>::New
(
    const meshLevel<Type,MeshType>& level,
    const boundary& b
)
{
    word boundaryConditionType;

    if (b.castable<patchBoundary>())
    {
        boundaryConditionType =
            level.field().found("boundaryConditions")
          ? word
            (
                level.field().subDict("boundaryConditions")
               .subDict(IOobject::member(b.name())).lookup("type")
            )
          : word("dummy");
    }
    else
    {
        boundaryConditionType = b.type();
    }

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
        cstrIter()(level, b)
    );
}

template<class Type, class MeshType>
inline labelVector boundaryCondition<Type,MeshType>::S(const label d) const
{
    return
        this->offset() == zeroXYZ
      ? zeroXYZ
      : this->fvMsh_.template S<MeshType>(l_,d,this->offset());
}

template<class Type, class MeshType>
inline labelVector boundaryCondition<Type,MeshType>::E(const label d) const
{
    return
        this->offset() == zeroXYZ
      ? zeroXYZ
      : this->fvMsh_.template E<MeshType>(l_,d,this->offset());
}

template<class Type, class MeshType>
inline labelVector boundaryCondition<Type,MeshType>::N(const label d) const
{
    return
        this->offset() == zeroXYZ
      ? zeroXYZ
      : this->fvMsh_.template N<MeshType>(l_,d,this->offset());
}

template<class Type, class MeshType>
const faceField<vector,MeshType>&
boundaryCondition<Type,MeshType>::faceCenters() const
{
    return fvMsh_.template metrics<MeshType>().faceCenters();
}

template<class Type, class MeshType>
const faceField<vector,MeshType>&
boundaryCondition<Type,MeshType>::faceNormals() const
{
    return fvMsh_.template metrics<MeshType>().faceNormals();
}

template<class Type, class MeshType>
const faceField<scalar,MeshType>&
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
const faceField<scalar,MeshType>&
boundaryCondition<Type,MeshType>::faceDeltas() const
{
    return fvMsh_.template metrics<MeshType>().faceDeltas();
}

// Level ghost eliminate

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::eliminateGhosts
(
    linearSystem<stencil,Type,MeshType>& sys
)
{
    if (this->eliminated())
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            eliminateGhosts(sys, d);
}

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::eliminateGhosts
(
    linearSystem<diagStencil,Type,MeshType>& sys
)
{
    if (this->eliminated())
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            eliminateGhosts(sys, d);
}

// Direction ghost eliminate

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::eliminateGhosts
(
    linearSystem<stencil,Type,MeshType>&,
    const label
)
{
    NotImplemented;
}

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::eliminateGhosts
(
    linearSystem<diagStencil,Type,MeshType>&,
    const label
)
{
    NotImplemented;
}

// Prepare

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::prepare()
{
    forAll(this->level_, d)
        prepare(d);
}

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::prepare(const label d)
{}

// Evaluate

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::evaluate()
{
    forAll(this->level_, d)
        evaluate(d);
}

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::evaluate(const label d)
{}

}

}

}
