#include "meshFieldSet.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::allocate
(
    const bool deep,
    const IOobject::readOption r,
    const IOobject::writeOption w
)
{
    listType::setSize(N);

    for (direction s = 0; s < N; s++)
    {
        listType::set
        (
            s,
            new meshField<Type,MeshType>
            (
                IOobject
                (
                    this->name() + "_" + Foam::name(s),
                    fvMsh_.time().path()/"0",
                    fvMsh_.db(),
                    r,
                    w,
                    this->registerObject()
                ),
                fvMsh_,
                deep
            )
        );
    }
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::transfer
(
    meshFieldSet<Type,MeshType,N>& set
)
{
    listType::transfer(set);
}

// Main constructors

template<class Type, class MeshType, direction N>
meshFieldSet<Type,MeshType,N>::meshFieldSet
(
    const word& name,
    const fvMesh& fvMsh,
    const IOobject::readOption r,
    const IOobject::writeOption w,
    const bool registerObject,
    const bool deep
)
:
    FastPtrList<meshField<Type,MeshType>>(),
    regIOobject
    (
        IOobject
        (
            name,
            fvMsh.time().path()/"0",
            fvMsh.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(fvMsh)
{
    this->allocate(deep, r, w);
}

template<class Type, class MeshType, direction N>
meshFieldSet<Type,MeshType,N>::meshFieldSet
(
    const IOobject& io,
    const fvMesh& fvMsh,
    const bool deep
)
:
    FastPtrList<meshField<Type,MeshType>>(),
    regIOobject
    (
        IOobject
        (
            io.name(),
            fvMsh.time().path()/"0",
            fvMsh.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            io.registerObject()
        )
    ),
    cachedRefCount(),
    fvMsh_(fvMsh)
{
    this->allocate(deep, io.readOpt(), io.writeOpt());
}

// Copy constructors

template<class Type, class MeshType, direction N>
meshFieldSet<Type,MeshType,N>::meshFieldSet
(
    const meshFieldSet<Type,MeshType,N>& set,
    const bool registerObject
)
:
    FastPtrList<meshField<Type,MeshType>>(set),
    regIOobject
    (
        IOobject
        (
            set.name(),
            set.fvMsh().time().path()/"0",
            set.fvMsh().db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(set.fvMsh_)
{}

template<class Type, class MeshType, direction N>
meshFieldSet<Type,MeshType,N>::meshFieldSet
(
    const word& name,
    const meshFieldSet<Type,MeshType,N>& set,
    const bool registerObject
)
:
    FastPtrList<meshField<Type,MeshType>>(set),
    regIOobject
    (
        IOobject
        (
            name,
            set.fvMsh().time().path()/"0",
            set.fvMsh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(set.fvMsh())
{
    rename(name);
}

template<class Type, class MeshType, direction N>
meshFieldSet<Type,MeshType,N>::meshFieldSet
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const bool registerObject
)
:
    FastPtrList<meshField<Type,MeshType>>
    (
        const_cast<meshFieldSet<Type,MeshType,N>&>(tSet()),
        tSet.isTmp() && tSet->unique()
    ),
    regIOobject
    (
        IOobject
        (
            tSet->name(),
            tSet->fvMsh_.time().path()/"0",
            tSet->fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(tSet->fvMsh_)
{
    if (tSet.isTmp())
        tSet.clear();
}

template<class Type, class MeshType, direction N>
meshFieldSet<Type,MeshType,N>::meshFieldSet
(
    const word& name,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const bool registerObject
)
:
    FastPtrList<meshField<Type,MeshType>>
    (
        const_cast<meshFieldSet<Type,MeshType,N>&>(tSet()),
        tSet.isTmp() && tSet->unique()
    ),
    regIOobject
    (
        IOobject
        (
            name,
            tSet->fvMsh_.time().path()/"0",
            tSet->fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(tSet->fvMsh_)
{
    rename(name);

    if (tSet.isTmp())
        tSet.clear();
}

template<class Type, class MeshType, direction N>
meshFieldSet<Type,MeshType,N>::~meshFieldSet()
{}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator=
(
    const meshFieldSet<Type,MeshType,N>& set
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) = set[s];
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator=
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    if (tSet.isTmp() && tSet->unique())
    {
        meshFieldSet<Type,MeshType,N>& set =
            const_cast<meshFieldSet<Type,MeshType,N>&>(tSet());

        transfer(set);

        tSet.clear();
    }
    else
    {
        *this = tSet();
    }
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator=
(
    const zero
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) = Zero;
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator=
(
    const Type& v
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) = v;
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator+=
(
    const Type& v
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) += v;
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator+=
(
    const meshFieldSet<Type,MeshType,N>& set
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) += set[s];
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator+=
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    *this += tSet();

    if (tSet.isTmp())
        tSet.clear();
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator-=
(
    const Type& v
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) -= v;
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator-=
(
    const meshFieldSet<Type,MeshType,N>& set
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) -= set[s];
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator-=
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    *this -= tSet();

    if (tSet.isTmp())
        tSet.clear();
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator*=
(
    const scalar s
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) *= s;
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator*=
(
    const meshFieldSet<scalar,MeshType,N>& set
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) *= set[s];
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator*=
(
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet
)
{
    *this *= tSet();

    if (tSet.isTmp())
        tSet.clear();
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator/=
(
    const scalar s
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) /= s;
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator/=
(
    const meshFieldSet<scalar,MeshType,N>& set
)
{
    for (direction s = 0; s < N; s++)
        listType::operator[](s) /= set[s];
}

template<class Type, class MeshType, direction N>
void meshFieldSet<Type,MeshType,N>::operator/=
(
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet
)
{
    *this /= tSet();

    if (tSet.isTmp())
        tSet.clear();
}

}

}

}
