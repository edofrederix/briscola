#include "boundaryCondition.H"

#include "colocated.H"
#include "staggered.H"

#include "stencil.H"
#include "diagStencil.H"

#include "meshField.H"

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
                if (bcs[bci].boundaryOffset() == bo)
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

template<class Type, class MeshType>
tmp<block<Type>> boundaryCondition<Type,MeshType>::internalValue
(
    const label l,
    const label d
)
{
    const labelVector bo(this->boundaryOffset());
    const meshDirection<Type,MeshType>& fld = mshField_[l][d];

    const labelVector S(fld.boundaryStart(bo));
    const labelVector E(fld.boundaryEnd(bo));

    return tmp<block<Type>>(new block<Type>(E-S, Zero));
}

template<class Type, class MeshType>
tmp<block<Type>> boundaryCondition<Type,MeshType>::boundarySources
(
    const label l,
    const label d
)
{
    const labelVector bo(this->boundaryOffset());
    const meshDirection<Type,MeshType>& fld = mshField_[l][d];

    const labelVector S(fld.boundaryStart(bo));
    const labelVector E(fld.boundaryEnd(bo));

    return tmp<block<Type>>(new block<Type>(E-S, Zero));
}

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::correctSystem
(
    linearSystem<stencil,Type,MeshType>& sys,
    const label l
)
{
    // Only face boundaries can manipulate a stencil

    if (this->boundaryOffsetDegree() == 1)
    {
        const labelVector bo(this->boundaryOffset());
        const label faceNum(faceNumber(bo));

        forAll(mshField_[l], d)
        {
            const meshDirection<Type,MeshType>& fld = this->mshField_[l][d];

            // Manipulate the linear system for eliminated or constrained
            // boundary conditions

            const bool shifted = fld.shifted(bo);
            const bool eliminated = this->eliminated(shifted);
            const bool constrained = this->constrained(shifted);

            if (eliminated || constrained)
            {
                meshDirection<stencil,MeshType>& Ad = sys.A()[l][d];
                meshDirection<Type,MeshType>& bd = sys.b()[l][d];

                const meshDirection<scalar,MeshType>& cv =
                    this->cellVolumes()[l][d];

                const labelVector S(fld.boundaryStart(bo));
                const labelVector E(fld.boundaryEnd(bo));

                labelVector ijk;

                if (eliminated)
                {
                    const stencil C(this->boundaryCoeff(l,d));
                    const block<Type> B(this->boundarySources(l,d));

                    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                    {
                        // Eliminate dependence on the ghost

                        Ad(ijk) += Ad(ijk)[faceNum+1]*C;
                        bd(ijk) -= Ad(ijk)[faceNum+1]*B(ijk-S);

                        Ad(ijk)[faceNum+1] = 0;
                    }
                }

                if (constrained)
                {
                    const block<Type> V(this->internalValue(l,d));

                    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                    {
                        // Constrain the internal value

                        Ad(ijk) = diagStencil(cv(ijk));
                        bd(ijk) = V(ijk-S)*cv(ijk);
                    }
                }
            }
        }
    }
}

template<class Type, class MeshType>
void boundaryCondition<Type,MeshType>::correctSystem
(
    linearSystem<diagStencil,Type,MeshType>& sys,
    const label l
)
{
    // Nothing to be done for a diagonal stencil
}

}

}

}
