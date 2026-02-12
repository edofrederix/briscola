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
boundaryConditionBaseType
boundaryCondition<Type,MeshType>::globalBaseType
(
    const meshLevel<Type,MeshType>& level,
    const labelVector bo
)
{
    // Global BCs only exist on rectilinear brick topologies

    const fvMesh& fvMsh = level.fvMsh();

    if (!fvMsh.msh().topology().rectilinear())
    {
        FatalErrorInFunction
            << "Brick topology is not rectlinear" << endl
            << abort(FatalError);

        return boundaryConditionBaseType::DUMMYBC;
    }
    else
    {
        const decompositionMap& map = level.lvl().decomp().map();

        const label facei = faceNumber(bo);
        const label dir = facei/2;
        const label i = -(facei%2);

        // Make sure that the boundary conditions are set

        if (!level.boundaryConditions().size())
            const_cast<meshLevel<Type,MeshType>&>(level)
                .addBoundaryConditions();

        const PtrList<boundaryCondition<Type,MeshType>>& bcs =
            level.boundaryConditions();

        // Block of processors at the requested face

        const labelBlock procs(map.slice(i,dir));

        labelList types(procs.size());

        forAllBlockLinear(procs, j)
        if (procs(j) == Pstream::myProcNo())
        {
            // Find type

            label type = -1;

            forAll(bcs, bci)
            {
                if (bcs[bci].offset() == bo)
                {
                    type = static_cast<label>(bcs[bci].baseType());
                    break;
                }
            }

            if (type == -1)
            {
                FatalErrorInFunction
                    << "Cannot find boundary condition type" << endl
                    << abort(FatalError);
            }

            if (procs(j) == Pstream::masterNo())
            {
                types[j] = type;
            }
            else
            {
                // Send type to master

                OPstream send
                (
                    Pstream::commsTypes::blocking,
                    Pstream::masterNo(),
                    0,
                    UPstream::msgType(),
                    level.lvl().comms()
                );

                send << type;
            }
        }

        // Receive

        if (Pstream::master())
        {
            forAllBlockLinear(procs, j)
            if (procs(j) != Pstream::masterNo())
            {
                IPstream recv
                (
                    Pstream::commsTypes::blocking,
                    procs(j),
                    0,
                    UPstream::msgType(),
                    level.lvl().comms()
                );

                recv >> types[j];
            }
        }

        // Check consistency

        label type;

        if (Pstream::master())
        {
            type = types[0];

            for (int j = 1; j < types.size(); j++)
            {
                if (types[j] != type)
                {
                    FatalErrorInFunction
                        << "Inconsistent boundary condition types for level "
                        << level.levelNum() << " at boundary offset "
                        << bo << endl << abort(FatalError);
                }
            }

            for (int proc = 0; proc < Pstream::nProcs(); proc++)
            if (proc != Pstream::masterNo())
            {
                OPstream send
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    0,
                    UPstream::msgType(),
                    level.lvl().comms()
                );

                send << type;
            }
        }
        else
        {
            IPstream recv
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo(),
                0,
                UPstream::msgType(),
                level.lvl().comms()
            );

            recv >> type;
        }

        return static_cast<boundaryConditionBaseType>(type);
    }
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
