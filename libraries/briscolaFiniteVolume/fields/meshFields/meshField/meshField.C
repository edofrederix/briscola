#include "meshField.H"
#include "fvMesh.H"

#include "domainBoundary.H"
#include "parallelBoundary.H"
#include "periodicBoundary.H"
#include "emptyBoundary.H"

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Private

template<class Type, class MeshType>
void meshField<Type,MeshType>::allocate(const bool deep)
{
    listType::setSize(deep ? fvMsh_.size() : 1);

    forAll(*this, l)
    {
        listType::set
        (
            l,
            new meshLevel<Type,MeshType>
            (
                *this,
                fvMsh_,
                l
            )
        );
    }
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::setFieldPointers()
{
    forAll(*this, l)
        listType::operator[](l).mshFieldPtr_ = this;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::transfer
(
    meshField<Type,MeshType>& field
)
{
    listType::transfer(field);

    oldTimePtr_ = field.oldTimePtr_;
    field.oldTimePtr_ = nullptr;

    setFieldPointers();
}

// Main constructor

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const word& name,
    const fvMesh& fvMsh,
    const IOobject::readOption r,
    const IOobject::writeOption w,
    const bool registerObject,
    const bool initBCs,
    const bool deep
)
:
    PtrList<meshLevel<Type,MeshType>>(),
    IOdictionary
    (
        IOobject
        (
            name,
            fvMsh.time().path()/"0",
            fvMsh.time(),
            r,
            w,
            registerObject
        )
    ),
    refCount(),
    fvMsh_(fvMsh),
    oldTimePtr_(nullptr),
    boundaryConditions_(),
    reScheme_()
{
    if (!fvMsh.structured() && MeshType::numberOfDirections > 1)
    {
        FatalErrorInFunction
            << "Cannot create a " << MeshType::typeName << " field on an "
            << "unstructured mesh."
            << endl << abort(FatalError);
    }

    if
    (
        (
            r == IOobject::MUST_READ
         || r == IOobject::MUST_READ_IF_MODIFIED
        )
     && !headerOk()
    )
    {
        FatalErrorInFunction
            << "Cannot read field from " << objectPath() << endl
            << exit(FatalError);
    }

    this->allocate(deep);

    if (initBCs)
        addBoundaryConditions();
}

// Copy constructors

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const meshField<Type,MeshType>& field,
    const bool registerObject,
    const bool copyBCs
)
:
    PtrList<meshLevel<Type,MeshType>>(field),
    IOdictionary
    (
        IOobject
        (
            field.name(),
            field.fvMsh().time().path()/"0",
            field.fvMsh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    refCount(),
    fvMsh_(field.fvMsh()),
    oldTimePtr_(nullptr),
    boundaryConditions_(),
    reScheme_()
{
    setFieldPointers();

    if (copyBCs)
    {
        PtrList<boundaryCondition<Type,MeshType>> list
        (
            field.boundaryConditions(),
            *this
        );

        boundaryConditions_.transfer(list);
    }
}

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const word& name,
    const meshField<Type,MeshType>& field,
    const bool registerObject,
    const bool copyBCs
)
:
    PtrList<meshLevel<Type,MeshType>>(field),
    IOdictionary
    (
        IOobject
        (
            name,
            field.fvMsh().time().path()/"0",
            field.fvMsh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    refCount(),
    fvMsh_(field.fvMsh()),
    oldTimePtr_(nullptr),
    boundaryConditions_(),
    reScheme_()
{
    setFieldPointers();

    if (copyBCs)
    {
        PtrList<boundaryCondition<Type,MeshType>> list
        (
            field.boundaryConditions(),
            *this
        );

        boundaryConditions_.transfer(list);
    }
}

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const tmp<meshField<Type,MeshType>>& tfield,
    const bool registerObject,
    const bool copyBCs
)
:
    PtrList<meshLevel<Type,MeshType>>
    (
        const_cast<meshField<Type,MeshType>&>(tfield()),
        tfield.isTmp()
    ),
    IOdictionary
    (
        IOobject
        (
            tfield->name(),
            tfield->fvMsh_.time().path()/"0",
            tfield->fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    refCount(),
    fvMsh_(tfield->fvMsh_),
    oldTimePtr_(),
    boundaryConditions_(),
    reScheme_()
{
    setFieldPointers();

    if (copyBCs)
    {
        PtrList<boundaryCondition<Type,MeshType>> list
        (
            tfield->boundaryConditions(),
            *this
        );

        boundaryConditions_.transfer(list);
    }

    tfield.clear();
}

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const word& name,
    const tmp<meshField<Type,MeshType>>& tfield,
    const bool registerObject,
    const bool copyBCs
)
:
    PtrList<meshLevel<Type,MeshType>>
    (
        const_cast<meshField<Type,MeshType>&>(tfield()),
        tfield.isTmp()
    ),
    IOdictionary
    (
        IOobject
        (
            name,
            tfield->fvMsh_.time().path()/"0",
            tfield->fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    refCount(),
    fvMsh_(tfield->fvMsh_),
    oldTimePtr_(),
    boundaryConditions_(),
    reScheme_()
{
    setFieldPointers();

    if (copyBCs)
    {
        PtrList<boundaryCondition<Type,MeshType>> list
        (
            tfield->boundaryConditions(),
            *this
        );

        boundaryConditions_.transfer(list);
    }

    tfield.clear();
}

// Destructor

template<class Type, class MeshType>
meshField<Type,MeshType>::~meshField()
{
    if (oldTimePtr_ != nullptr)
    {
         delete oldTimePtr_;
    }
}

// Public

template<class Type, class MeshType>
void meshField<Type,MeshType>::addBoundaryConditions()
{
    if (boundaryConditions_.size() == 0)
    {
        // First add domain boundaries

        forAll(fvMsh_.boundaries(), bi)
        {
            const boundary& b = fvMsh_.boundaries()[bi];

            if (b.castable<domainBoundary>())
            {
                boundaryConditions_.append
                (
                    boundaryCondition<Type,MeshType>::NewDomain
                    (
                        *this,
                        b
                    )
                );
            }
        }

        // Add parallel and periodic boundaries. First faces, then edges and
        // finally vertices.

        for(label order = 1; order <= 3; order++)
        forAll(fvMsh_.boundaries(), bi)
        {
            const boundary& b = fvMsh_.boundaries()[bi];
            const labelVector bo(b.offset());

            if (cmptSum(cmptMag(bo)) == order)
            {
                if (b.castable<periodicBoundary>())
                {
                    boundaryConditions_.append
                    (
                        boundaryCondition<Type,MeshType>::NewPeriodic
                        (
                            *this,
                            b
                        )
                    );
                }
                else if (b.castable<parallelBoundary>())
                {
                    boundaryConditions_.append
                    (
                        boundaryCondition<Type,MeshType>::NewParallel
                        (
                            *this,
                            b
                        )
                    );
                }
            }
        }

        // Finally add empties

        forAll(fvMsh_.boundaries(), bi)
        {
            const boundary& b = fvMsh_.boundaries()[bi];

            if (b.castable<emptyBoundary>())
            {
                boundaryConditions_.append
                (
                    boundaryCondition<Type,MeshType>::NewEmpty
                    (
                        *this,
                        b
                    )
                );
            }
        }
    }
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::correctBoundaryConditions()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctBoundaryConditions();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::correctParallelBoundaryConditions()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctParallelBoundaryConditions();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::correctPeriodicBoundaryConditions()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctPeriodicBoundaryConditions();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::correctNonEliminatedBoundaryConditions()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctNonEliminatedBoundaryConditions();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::correctEliminatedBoundaryConditions()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctEliminatedBoundaryConditions();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::setOldTime()
{
    if (oldTimePtr_ == nullptr)
    {
        oldTimePtr_ = new meshField<Type,MeshType>
        (
            IOobject::groupName(name(), "oldTime"),
            fvMsh_
        );
    }

    *oldTimePtr_ = *this;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::makeDeep()
{
    if (this->shallow())
    {
        listType::setSize(fvMsh_.size());

        for (int l = 1; l < fvMsh_.size(); l++)
        {
            listType::set
            (
                l,
                new meshLevel<Type,MeshType>
                (
                    *this,
                    fvMsh_,
                    l
                )
            );
        }
    }
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::makeShallow()
{
    if (this->deep())
        listType::setSize(1);
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::setRestrictionScheme(const word scheme)
{
    reScheme_.reset
    (
        restrictionScheme<Type,MeshType>::New(scheme, fvMsh_).ptr()
    );
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::restrict()
{
    if (reScheme_.empty())
        this->setRestrictionScheme
        (
            restrictionScheme<Type,MeshType>::defaultScheme
        );

    makeDeep();

    reScheme_->restrict(*this);
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator=(const meshField<Type,MeshType>& F)
{
    this->make(F.deep());
    forAll(*this, l)
        listType::operator[](l) = F[l];
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator=
(
    const tmp<meshField<Type,MeshType>>& tF
)
{
    if (tF.isTmp())
    {
        meshField<Type,MeshType>& F =
            const_cast<meshField<Type,MeshType>&>(tF());

        transfer(F);
    }
    else
    {
        this->make(tF->deep());
        *this = tF();
    }

    tF.clear();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator=(const Type& v)
{
    forAll(*this, l)
        listType::operator[](l) = v;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator=(const List<Type>& v)
{
    forAll(*this, l)
        listType::operator[](l) = v;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator=(const zero)
{
    forAll(*this, l)
        listType::operator[](l) = Zero;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator+=(const meshField<Type,MeshType>& F)
{
    checkMemberOperatorArgDepth(F);
    forAll(*this, l)
        listType::operator[](l) += F[l];
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator+=
(
    const tmp<meshField<Type,MeshType>>& tF
)
{
    checkMemberOperatorArgDepth(tF());
    *this += tF();
    tF.clear();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator-=(const meshField<Type,MeshType>& F)
{
    checkMemberOperatorArgDepth(F);
    forAll(*this, l)
        listType::operator[](l) -= F[l];
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator-=
(
    const tmp<meshField<Type,MeshType>>& tF
)
{
    checkMemberOperatorArgDepth(tF());
    *this -= tF();
    tF.clear();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator*=(const meshField<scalar,MeshType>& F)
{
    checkMemberOperatorArgDepth(F);
    forAll(*this, l)
        listType::operator[](l) *= F[l];
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator*=
(
    const tmp<meshField<scalar,MeshType>>& tF
)
{
    checkMemberOperatorArgDepth(tF());
    *this *= tF();
    tF.clear();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator/=(const meshField<scalar,MeshType>& F)
{
    checkMemberOperatorArgDepth(F);
    forAll(*this, l)
        listType::operator[](l) /= F[l];
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator/=
(
    const tmp<meshField<scalar,MeshType>>& tF
)
{
    checkMemberOperatorArgDepth(tF());
    *this /= tF();
    tF.clear();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator+=(const Type& v)
{
    forAll(*this, l)
        listType::operator[](l) += v;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator+=(const List<Type>& v)
{
    forAll(*this, l)
        listType::operator[](l) += v;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator-=(const Type& v)
{
    forAll(*this, l)
        listType::operator[](l) -= v;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator-=(const List<Type>& v)
{
    forAll(*this, l)
        listType::operator[](l) -= v;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator*=(const scalar& v)
{
    forAll(*this, l)
        listType::operator[](l) *= v;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator*=(const scalarList& v)
{
    forAll(*this, l)
        listType::operator[](l) *= v;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator/=(const scalar& v)
{
    forAll(*this, l)
        listType::operator[](l) /= v;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator/=(const scalarList& v)
{
    forAll(*this, l)
        listType::operator[](l) /= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator=
(
    const meshField<Type2,MeshType>& F
)
{
    this->make(F.deep());
    forAll(*this, l)
        listType::operator[](l) = F[l];
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator=
(
    const tmp<meshField<Type2,MeshType>>& tF
)
{
    this->make(tF->deep());
    *this = tF();
    tF.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator+=
(
    const meshField<Type2,MeshType>& F
)
{
    checkMemberOperatorArgDepth(F);
    forAll(*this, l)
        listType::operator[](l) += F[l];
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator+=
(
    const tmp<meshField<Type2,MeshType>>& tF
)
{
    checkMemberOperatorArgDepth(tF());
    *this += tF();
    tF.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator-=
(
    const meshField<Type2,MeshType>& F
)
{
    checkMemberOperatorArgDepth(F);
    forAll(*this, l)
        listType::operator[](l) -= F[l];
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator-=
(
    const tmp<meshField<Type2,MeshType>>& tF
)
{
    checkMemberOperatorArgDepth(tF());
    *this -= tF();
    tF.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator*=
(
    const meshField<Type2,MeshType>& F
)
{
    checkMemberOperatorArgDepth(F);
    forAll(*this, l)
        listType::operator[](l) *= F[l];
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator*=
(
    const tmp<meshField<Type2,MeshType>>& tF
)
{
    checkMemberOperatorArgDepth(tF());
    *this *= tF();
    tF.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator/=
(
    const meshField<Type2,MeshType>& F
)
{
    checkMemberOperatorArgDepth(F);
    forAll(*this, l)
        listType::operator[](l) /= F[l];
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator/=
(
    const tmp<meshField<Type2,MeshType>>& tF
)
{
    checkMemberOperatorArgDepth(tF());
    *this /= tF();
    tF.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator=(const Type2& v)
{
    forAll(*this, l)
        listType::operator[](l) = v;
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator=(const List<Type2>& v)
{
    forAll(*this, l)
        listType::operator[](l) = v;
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator+=(const Type2& v)
{
    forAll(*this, l)
        listType::operator[](l) += v;
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator+=(const List<Type2>& v)
{
    forAll(*this, l)
        listType::operator[](l) += v;
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator-=(const Type2& v)
{
    forAll(*this, l)
        listType::operator[](l) -= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator-=(const List<Type2>& v)
{
    forAll(*this, l)
        listType::operator[](l) -= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator*=(const Type2& v)
{
    forAll(*this, l)
        listType::operator[](l) *= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator*=(const List<Type2>& v)
{
    forAll(*this, l)
        listType::operator[](l) *= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator/=(const Type2& v)
{
    forAll(*this, l)
        listType::operator[](l) /= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator/=(const List<Type2>& v)
{
    forAll(*this, l)
        listType::operator[](l) /= v;
}

}

}

}

#include "meshFieldFunctions.C"
#include "meshFieldStencilFunctions.C"
