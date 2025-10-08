#include "boundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshField.H"
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
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    refCount(),
    fvMsh_(mshField.fvMsh()),
    mshField_(const_cast<meshField<Type,MeshType>&>(mshField)),
    b_(b),
    dict_
    (
        mshField_.found("boundaryConditions")
     && mshField_.subDict("boundaryConditions").found(IOobject::member(b_.name()))
      ? mshField_.subDict("boundaryConditions").subDict(IOobject::member(b_.name()))
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
    b_(bc.b_),
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
    b_(bc.b_),
    dict_(bc.dict_)
{}

template<class Type, class MeshType>
boundaryCondition<Type,MeshType>::~boundaryCondition()
{}

template<class Type, class MeshType>
autoPtr<boundaryCondition<Type,MeshType>> boundaryCondition<Type,MeshType>::New
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
{
    word boundaryConditionType;

    if (b.castable<patchBoundary>())
    {
        boundaryConditionType =
            mshField.found("boundaryConditions")
          ? word
            (
                mshField.subDict("boundaryConditions")
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
        cstrIter()(mshField, b)
    );
}

template<class Type, class MeshType>
boundaryConditionBaseType
boundaryCondition<Type,MeshType>::globalBaseType
(
    const meshField<Type,MeshType>& field,
    const labelVector bo
)
{
    // Global BCs only exist on rectilinear brick topologies

    const fvMesh& fvMsh = field.fvMsh();

    if (!fvMsh.msh().topology().rectilinear())
    {
        FatalErrorInFunction
            << "Brick topology is not rectlinear" << endl
            << abort(FatalError);

        return boundaryConditionBaseType::DUMMYBC;
    }
    else
    {
        const decompositionMap& map = fvMsh.msh().decomp().map();

        const label facei = faceNumber(bo);
        const label dir = facei/2;
        const label i = -(facei%2);

        // Make sure that the boundary conditions are set

        if (!field.boundaryConditions().size())
            const_cast<meshField<Type,MeshType>&>(field)
                .addBoundaryConditions();

        const PtrList<boundaryCondition<Type,MeshType>>& bcs =
            field.boundaryConditions();

        // Block of processors at the requested face

        const labelBlock procs(map.slice(i,dir));

        labelList types(procs.size());

        forAllBlockLinear(procs, i)
        if (procs(i) == Pstream::myProcNo())
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

            if (procs(i) == Pstream::masterNo())
            {
                types[i] = type;
            }
            else
            {
                // Send type to master

                OPstream send
                (
                    Pstream::commsTypes::blocking,
                    Pstream::masterNo()
                );

                send << type;
            }
        }

        // Receive

        if (Pstream::master())
        {
            forAllBlockLinear(procs, i)
            if (procs(i) != Pstream::masterNo())
            {
                IPstream recv
                (
                    Pstream::commsTypes::blocking,
                    procs(i)
                );

                recv >> types[i];
            }
        }

        // Check consistency

        label type;

        if (Pstream::master())
        {
            type = types[0];

            for (int i = 1; i < types.size(); i++)
            {
                if (types[i] != type)
                {
                    FatalErrorInFunction
                        << "Inconsistent boundary condition types for field "
                        << field.name() << " at boundary offset " << bo << endl
                        << abort(FatalError);
                }
            }

            for (int proc = 0; proc < Pstream::nProcs(); proc++)
            if (proc != Pstream::masterNo())
            {
                OPstream send
                (
                    Pstream::commsTypes::blocking,
                    proc
                );

                send << type;
            }
        }
        else
        {
            IPstream recv
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            recv >> type;
        }

        return static_cast<boundaryConditionBaseType>(type);
    }
}

template<class Type, class MeshType>
inline labelVector boundaryCondition<Type,MeshType>::S
(
    const label l,
    const label d
) const
{
    return
        this->offset() == zeroXYZ
      ? zeroXYZ
      : this->fvMsh_.template S<MeshType>(l,d,this->offset());
}

template<class Type, class MeshType>
inline labelVector boundaryCondition<Type,MeshType>::E
(
    const label l,
    const label d
) const
{
    return
        this->offset() == zeroXYZ
      ? zeroXYZ
      : this->fvMsh_.template E<MeshType>(l,d,this->offset());
}

template<class Type, class MeshType>
inline labelVector boundaryCondition<Type,MeshType>::N
(
    const label l,
    const label d
) const
{
    return
        this->offset() == zeroXYZ
      ? zeroXYZ
      : this->fvMsh_.template N<MeshType>(l,d,this->offset());
}

template<class Type, class MeshType>
const meshField<faceVector,MeshType>&
boundaryCondition<Type,MeshType>::faceCenters() const
{
    return fvMsh_.template metrics<MeshType>().faceCenters();
}

template<class Type, class MeshType>
const meshField<faceVector,MeshType>&
boundaryCondition<Type,MeshType>::faceNormals() const
{
    return fvMsh_.template metrics<MeshType>().faceNormals();
}

template<class Type, class MeshType>
const meshField<faceScalar,MeshType>&
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
const meshField<faceScalar,MeshType>&
boundaryCondition<Type,MeshType>::faceDeltas() const
{
    return fvMsh_.template metrics<MeshType>().faceDeltas();
}

// Level ghost eliminate

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::eliminateGhosts
(
    linearSystem<stencil,Type,MeshType>& sys,
    const label l
)
{
    if (this->eliminated())
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            eliminateGhosts(sys, l, d);
}

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::eliminateGhosts
(
    linearSystem<diagStencil,Type,MeshType>& sys,
    const label l
)
{
    if (this->eliminated())
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            eliminateGhosts(sys, l, d);
}

// Direction ghost eliminate

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::eliminateGhosts
(
    linearSystem<stencil,Type,MeshType>&,
    const label,
    const label
)
{
    NotImplemented;
}

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::eliminateGhosts
(
    linearSystem<diagStencil,Type,MeshType>&,
    const label,
    const label
)
{
    NotImplemented;
}

// Prepare

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::prepare(const label l)
{
    forAll(this->mshField_[l], d)
        prepare(l,d);
}

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::prepare(const label l, const label d)
{}

// Evaluate

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::evaluate(const label l)
{
    forAll(this->mshField_[l], d)
        evaluate(l,d);
}

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::evaluate(const label l, const label d)
{}

}

}

}
