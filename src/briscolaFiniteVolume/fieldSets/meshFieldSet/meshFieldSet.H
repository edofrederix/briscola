#ifndef meshFieldSet_H
#define meshFieldSet_H

#include "meshFieldSetFwd.H"
#include "meshFields.H"
#include "linearSystemFwd.H"
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
    public IOdictionary,
    public cachedRefCount
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Restriction scheme name
    word restrictionScheme_;

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
        const fvMesh& fvMsh
    )
    {
        tmp<meshFieldSet<Type,MeshType,N>> tSet =
            fvMsh.template
            NewObject<meshFieldSet<Type,MeshType,N>>(name);

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
        IOdictionary::rename(newName);

        forAll(*this, s)
            listType::operator[](s).rename(newName + "_" + Foam::name(s));
    }

    // Set restriction scheme
    void setRestrictionScheme(const word scheme)
    {
        restrictionScheme_ = scheme;
    }

    // Return restriction scheme name
    const word& restrictionScheme() const { return restrictionScheme_; }

    // Number of fields
    inline label size() const { return listType::size(); }

    // Access to the fields
    inline meshField<Type,MeshType>& operator[](label s)
    {
        return listType::operator[](s);
    }

    // Const access to the fields
    inline const meshField<Type,MeshType>& operator[](label s) const
    {
        return listType::operator[](s);
    }

    // Set maximum/minimum of field set and given value
    inline void max(const Type&);
    inline void min(const Type&);

    // Member Operators

    inline void operator=(const meshFieldSet<Type,MeshType,N>&);
    inline void operator=(const tmp<meshFieldSet<Type,MeshType,N>>&);
    inline void operator=(const Type&);
    inline void operator=(const zero);

    inline void operator+=(const Type&);
    inline void operator+=(const meshFieldSet<Type,MeshType,N>&);
    inline void operator+=(const tmp<meshFieldSet<Type,MeshType,N>>&);

    inline void operator-=(const Type&);
    inline void operator-=(const meshFieldSet<Type,MeshType,N>&);
    inline void operator-=(const tmp<meshFieldSet<Type,MeshType,N>>&);

    inline void operator*=(const scalar);
    inline void operator*=(const meshFieldSet<scalar,MeshType,N>&);
    inline void operator*=(const tmp<meshFieldSet<scalar,MeshType,N>>&);

    inline void operator/=(const scalar);
    inline void operator/=(const meshFieldSet<scalar,MeshType,N>&);
    inline void operator/=(const tmp<meshFieldSet<scalar,MeshType,N>>&);

    // Befriend linear system
    template<class SType, class, class>
    friend class linearSystem;

    // Befriend mesh level
    template<class, class>
    friend class meshField;
};

}

}

}

#include "meshFieldSet.C"
#include "meshFieldSetI.H"

#endif
