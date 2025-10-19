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
        // Fields are also registered if the set is registered. This is needed
        // for IO, since IO occurs through the fields, not the set.

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
    IOdictionary
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
    IOdictionary
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.db(),
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
    IOdictionary
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
    IOdictionary
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
    IOdictionary
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
    IOdictionary
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

}

}

}
