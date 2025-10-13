#include "meshField.H"
#include "fvMesh.H"

#include "patchBoundary.H"
#include "parallelBoundary.H"
#include "periodicBoundary.H"
#include "emptyBoundary.H"

#include "restrictionScheme.H"
#include "boundaryExchange.H"

#include "immersedBoundaryCondition.H"

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
            new meshLevel<Type,MeshType>(*this, l)
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
            fvMsh.db(),
            r,
            w,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(fvMsh),
    oldTimePtr_(nullptr),
    boundaryConditions_(),
    reSchemePtr_()
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
}

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const IOobject& io,
    const fvMesh& fvMsh,
    const bool deep
)
:
    PtrList<meshLevel<Type,MeshType>>(),
    IOdictionary(io),
    cachedRefCount(),
    fvMsh_(fvMsh),
    oldTimePtr_(nullptr),
    boundaryConditions_(),
    reSchemePtr_()
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
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     && !headerOk()
    )
    {
        FatalErrorInFunction
            << "Cannot read field from " << objectPath() << endl
            << exit(FatalError);
    }

    this->allocate(deep);
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
    cachedRefCount(),
    fvMsh_(field.fvMsh()),
    oldTimePtr_(nullptr),
    boundaryConditions_(),
    reSchemePtr_()
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
    cachedRefCount(),
    fvMsh_(field.fvMsh()),
    oldTimePtr_(nullptr),
    boundaryConditions_(),
    reSchemePtr_()
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
    const tmp<meshField<Type,MeshType>>& tField,
    const bool registerObject,
    const bool copyBCs
)
:
    PtrList<meshLevel<Type,MeshType>>
    (
        const_cast<meshField<Type,MeshType>&>(tField()),
        tField.isTmp() && tField->unique()
    ),
    IOdictionary
    (
        IOobject
        (
            tField->name(),
            tField->fvMsh_.time().path()/"0",
            tField->fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(tField->fvMsh_),
    oldTimePtr_(nullptr),
    boundaryConditions_(),
    reSchemePtr_()
{
    setFieldPointers();

    if (copyBCs)
    {
        PtrList<boundaryCondition<Type,MeshType>> list
        (
            tField->boundaryConditions(),
            *this
        );

        boundaryConditions_.transfer(list);
    }

    if (tField.isTmp())
        tField.clear();
}

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const word& name,
    const tmp<meshField<Type,MeshType>>& tField,
    const bool registerObject,
    const bool copyBCs
)
:
    PtrList<meshLevel<Type,MeshType>>
    (
        const_cast<meshField<Type,MeshType>&>(tField()),
        tField.isTmp() && tField->unique()
    ),
    IOdictionary
    (
        IOobject
        (
            name,
            tField->fvMsh_.time().path()/"0",
            tField->fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(tField->fvMsh_),
    oldTimePtr_(nullptr),
    boundaryConditions_(),
    reSchemePtr_()
{
    setFieldPointers();

    if (copyBCs)
    {
        PtrList<boundaryCondition<Type,MeshType>> list
        (
            tField->boundaryConditions(),
            *this
        );

        boundaryConditions_.transfer(list);
    }

    if (tField.isTmp())
        tField.clear();
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
        boundaryConditions_.setSize(fvMsh_.boundaries().size());

        forAll(fvMsh_.boundaries(), i)
        {
            boundaryConditions_.set
            (
                i,
                boundaryCondition<Type,MeshType>::New
                (
                    *this,
                    fvMsh_.boundaries()[i]
                )
            );
        }
    }

    #ifdef BOUNDARYEXCHANGE

    if (bExchangePtr_.empty())
        bExchangePtr_.reset(new boundaryExchange<Type,MeshType>(*this));

    #endif
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::addImmersedBoundaryConditions()
{
    if
    (
        !immersedBoundaryConditions_.size()
     && fvMsh_.immersedBoundaries<MeshType>().size()
    )
    {
        immersedBoundaryConditions_.setSize
        (
            fvMsh_.immersedBoundaries<MeshType>().size()
        );

        forAll(fvMsh_.immersedBoundaries<MeshType>(), i)
        {
            immersedBoundaryConditions_.set
            (
                i,
                immersedBoundaryCondition<Type,MeshType>::New
                (
                    *this,
                    fvMsh_.immersedBoundaries<MeshType>()[i]
                )
            );
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
void meshField<Type,MeshType>::correctPatchBoundaryConditions()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctPatchBoundaryConditions();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::correctEmptyBoundaryConditions()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctEmptyBoundaryConditions();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::correctCommsBoundaryConditions()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctCommsBoundaryConditions();
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
void meshField<Type,MeshType>::correctImmersedBoundaryConditions()
{
    addImmersedBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctImmersedBoundaryConditions();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::setOldTime()
{
    if (oldTimePtr_ == nullptr)
    {
        oldTimePtr_ = new meshField<Type,MeshType>
        (
            IOobject::groupName(name(), "oldTime"),
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            this->registerObject()
        );
    }

    if (oldTimePtr_->oldTimePtr_ != nullptr)
    {
        oldTimePtr_->setOldTime();
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
                new meshLevel<Type,MeshType>(*this, l)
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
    reSchemePtr_.reset
    (
        restrictionScheme<Type,MeshType>::New(fvMsh_, scheme).ptr()
    );
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::restrict()
{
    if (reSchemePtr_.empty())
        this->setRestrictionScheme
        (
            restrictionScheme<Type,MeshType>::defaultScheme
        );

    makeDeep();

    reSchemePtr_->restrict(*this);
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
    if (tF.isTmp() && tF->unique())
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

    if (tF.isTmp())
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
    if (F.shallow())
        makeShallow();

    forAll(*this, l)
        listType::operator[](l) += F[l];
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator+=
(
    const tmp<meshField<Type,MeshType>>& tF
)
{
    if (tF->shallow())
        makeShallow();

    *this += tF();
    if (tF.isTmp())
        tF.clear();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator-=(const meshField<Type,MeshType>& F)
{
    if (F.shallow())
        makeShallow();

    forAll(*this, l)
        listType::operator[](l) -= F[l];
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator-=
(
    const tmp<meshField<Type,MeshType>>& tF
)
{
    if (tF->shallow())
        makeShallow();

    *this -= tF();
    if (tF.isTmp())
        tF.clear();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator*=(const meshField<scalar,MeshType>& F)
{
    if (F.shallow())
        makeShallow();

    forAll(*this, l)
        listType::operator[](l) *= F[l];
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator*=
(
    const tmp<meshField<scalar,MeshType>>& tF
)
{
    if (tF->shallow())
        makeShallow();

    *this *= tF();
    if (tF.isTmp())
        tF.clear();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator/=(const meshField<scalar,MeshType>& F)
{
    if (F.shallow())
        makeShallow();

    forAll(*this, l)
        listType::operator[](l) /= F[l];
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::operator/=
(
    const tmp<meshField<scalar,MeshType>>& tF
)
{
    if (tF->shallow())
        makeShallow();

    *this /= tF();
    if (tF.isTmp())
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
    if (tF.isTmp())
        tF.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator+=
(
    const meshField<Type2,MeshType>& F
)
{
    if (F.shallow())
        makeShallow();

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
    if (tF->shallow())
        makeShallow();

    *this += tF();
    if (tF.isTmp())
        tF.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator-=
(
    const meshField<Type2,MeshType>& F
)
{
    if (F.shallow())
        makeShallow();

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
    if (tF->shallow())
        makeShallow();

    *this -= tF();
    if (tF.isTmp())
        tF.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator*=
(
    const meshField<Type2,MeshType>& F
)
{
    if (F.shallow())
        makeShallow();

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
    if (tF->shallow())
        makeShallow();

    *this *= tF();
    if (tF.isTmp())
        tF.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshField<Type,MeshType>::operator/=
(
    const meshField<Type2,MeshType>& F
)
{
    if (F.shallow())
        makeShallow();

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
    if (tF->shallow())
        makeShallow();

    *this /= tF();
    if (tF.isTmp())
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
