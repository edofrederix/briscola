#include "meshField.H"
#include "fvMesh.H"

#include "patchBoundary.H"
#include "parallelBoundary.H"
#include "periodicBoundary.H"
#include "emptyBoundary.H"

#include "restrictionScheme.H"

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
        listType::operator[](l).fieldPtr_ = this;
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
    FastPtrList<meshLevel<Type,MeshType>>(),
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
    oldTimePtr_(nullptr)
{
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
    FastPtrList<meshLevel<Type,MeshType>>(),
    IOdictionary(io),
    cachedRefCount(),
    fvMsh_(fvMsh),
    oldTimePtr_(nullptr)
{
    this->allocate(deep);
}

// Copy constructors

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const meshField<Type,MeshType>& field,
    const bool registerObject
)
:
    FastPtrList<meshLevel<Type,MeshType>>(field),
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
    oldTimePtr_(nullptr)
{
    setFieldPointers();
}

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const word& name,
    const meshField<Type,MeshType>& field,
    const bool registerObject
)
:
    FastPtrList<meshLevel<Type,MeshType>>(field),
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
    oldTimePtr_(nullptr)
{
    setFieldPointers();
}

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const tmp<meshField<Type,MeshType>>& tField,
    const bool registerObject
)
:
    FastPtrList<meshLevel<Type,MeshType>>
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
    oldTimePtr_(nullptr)
{
    setFieldPointers();

    if (tField.isTmp())
        tField.clear();
}

template<class Type, class MeshType>
meshField<Type,MeshType>::meshField
(
    const word& name,
    const tmp<meshField<Type,MeshType>>& tField,
    const bool registerObject
)
:
    FastPtrList<meshLevel<Type,MeshType>>
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
    oldTimePtr_(nullptr)
{
    setFieldPointers();

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
    forAll(*this, i)
        this->operator[](i).addBoundaryConditions();
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
void meshField<Type,MeshType>::transferBoundaryConditions
(
    const meshField<Type,MeshType>& field
)
{
    forAll(*this, i)
        this->operator[](i).transferBoundaryConditions(field[i]);
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::transferImmersedBoundaryConditions
(
    const meshField<Type,MeshType>& F
)
{
    PtrList<immersedBoundaryCondition<Type,MeshType>> list
    (
        F.immersedBoundaryConditions(),
        *this
    );

    immersedBoundaryConditions_.clear();
    immersedBoundaryConditions_.transfer(list);
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::correctBoundaryConditions()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctBoundaryConditions();
}

template<class Type, class MeshType>
template<class Selector>
void meshField<Type,MeshType>::prepare()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).template prepare<Selector>();
}

template<class Type, class MeshType>
template<class Selector>
void meshField<Type,MeshType>::evaluate()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).template evaluate<Selector>();
}

template<class Type, class MeshType>
template<class Selector>
void meshField<Type,MeshType>::correct()
{
    addBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).template correct<Selector>();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::correctImmersedBoundaryConditions()
{
    addImmersedBoundaryConditions();

    forAll(*this, l)
        listType::operator[](l).correctImmersedBoundaryConditions();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::correctAggData()
{
    forAll(*this, l)
        listType::operator[](l).correctAggData();
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
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
meshField<Type,MeshType>::component
(
    const label dir
) const
{
    tmp<meshField<cmptType,MeshType>> tF =
        meshField<cmptType,MeshType>::New
        (
            this->name() + ".component(" + Foam::name(dir) + ")",
            this->fvMsh_
        );

    tF.ref().make(deep());

    forAll(*this, l)
        tF.ref()[l] = listType::operator[](l).component(dir);

    return tF;
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::replace
(
    const label dir,
    const List<cmptType>& values
)
{
    forAll(*this, l)
        listType::operator[](l).replace(dir,values);
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::replace
(
    const label dir,
    const meshField<cmptType,MeshType>& F
)
{
    forAll(*this, l)
        listType::operator[](l).replace(dir,F[l]);
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::replace
(
    const label dir,
    const tmp<meshField<cmptType,MeshType>>& tF
)
{
    this->replace(dir,tF());

    if (tF.isTmp())
        tF.clear();
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::max(const Type& v)
{
    forAll(*this, l)
        listType::operator[](l).max(v);
}

template<class Type, class MeshType>
void meshField<Type,MeshType>::min(const Type& v)
{
    forAll(*this, l)
        listType::operator[](l).min(v);
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
