#ifndef meshFieldSet_H
#define meshFieldSet_H

#include "meshFieldSetFwd.H"
#include "meshFields.H"
#include "cachedRefCount.H"
#include "direction.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType, direction N>
class meshFieldSet
:
    public FastPtrList<meshField<Type,MeshType>>,
    public regIOobject,
    public cachedRefCount
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Allocate data and specify field depth and read/write options
    void allocate
    (
        const bool,
        const IOobject::readOption,
        const IOobject::writeOption
    );

    // Transfer data from other set
    void transfer(meshFieldSet<Type,MeshType,N>&);

public:

    // Number of components in the set
    static const direction nSetComponents = N;

    // Data type
    typedef Type dataType;

    // Mesh type
    typedef MeshType meshType;

    // Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    // List type
    typedef FastPtrList<meshField<Type,MeshType>> listType;

    // Main constructor: construct from name and mesh with optional read and
    // write options, object registration and depth

    meshFieldSet
    (
        const word& name,
        const fvMesh& fvMsh,
        const IOobject::readOption r = NO_READ,
        const IOobject::writeOption w = NO_WRITE,
        const bool registerObjects = false,
        const bool deep = false
    );

    meshFieldSet
    (
        const IOobject& io,
        const fvMesh& fvMsh,
        const bool deep = false
    );

    // Copy constructors: construct from set or tmp of set

    meshFieldSet
    (
        const meshFieldSet<Type,MeshType,N>& set,
        const bool registerObject = false
    );

    meshFieldSet
    (
        const word& name,
        const meshFieldSet<Type,MeshType,N>& set,
        const bool registerObject = false
    );

    meshFieldSet
    (
        const tmp<meshFieldSet<Type,MeshType,N>>& set,
        const bool registerObject = false
    );

    meshFieldSet
    (
        const word& name,
        const tmp<meshFieldSet<Type,MeshType,N>>& set,
        const bool registerObject = false
    );

    // Temporary set constructor
    static tmp<meshFieldSet<Type,MeshType,N>> New
    (
        const word& name,
        const fvMesh& fvMsh,
        const bool cached = false
    )
    {
        tmp<meshFieldSet<Type,MeshType,N>> tSet =
            fvMsh.template
            NewObject<meshFieldSet<Type,MeshType,N>>(name, cached);

        tSet.ref().rename(name);

        return tSet;
    }

    // Clone
    inline autoPtr<meshFieldSet<Type,MeshType,N>> clone() const
    {
        return autoPtr<meshFieldSet<Type,MeshType,N>>
        (
            new meshFieldSet<Type,MeshType,N>(*this)
        );
    }

    // Destructor
    virtual ~meshFieldSet();

    // Return the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Write data (not supported)
    bool writeData(Ostream& os) const { return os; }

    // When renaming, make sure to also rename the fields in the set
    virtual void rename(const word& newName)
    {
        IOobject::rename(newName);

        for (direction s = 0; s < N; s++)
            listType::operator[](s).rename(newName + "_" + Foam::name(s));
    }

    // Member Operators

    void operator=(const meshFieldSet<Type,MeshType,N>&);
    void operator=(const tmp<meshFieldSet<Type,MeshType,N>>&);
    void operator=(const Type&);
    void operator=(const zero);

    void operator+=(const Type&);
    void operator+=(const meshFieldSet<Type,MeshType,N>&);
    void operator+=(const tmp<meshFieldSet<Type,MeshType,N>>&);

    void operator-=(const Type&);
    void operator-=(const meshFieldSet<Type,MeshType,N>&);
    void operator-=(const tmp<meshFieldSet<Type,MeshType,N>>&);

    void operator*=(const scalar);
    void operator*=(const meshFieldSet<scalar,MeshType,N>&);
    void operator*=(const tmp<meshFieldSet<scalar,MeshType,N>>&);

    void operator/=(const scalar);
    void operator/=(const meshFieldSet<scalar,MeshType,N>&);
    void operator/=(const tmp<meshFieldSet<scalar,MeshType,N>>&);
};

}

}

}

#include "meshFieldSetI.H"

#ifdef NoRepository
    #include "meshFieldSet.C"
#endif

#endif
