#ifndef meshField_H
#define meshField_H

#include "block.H"
#include "IOdictionary.H"
#include "autoPtr.H"

#include "meshLevel.H"
#include "boundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh;

template<class SType, class Type, class MeshType>
class linearSystem;

template<class Type, class MeshType>
class meshField
:
    public PtrList<meshLevel<Type,MeshType>>,
    public IOdictionary,
    public refCount
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Old time field pointer, explicitly set by setOldTime and never changed by
    // any operation
    meshField<Type,MeshType>* oldTimePtr_;

    // Boundary conditions
    PtrList<boundaryCondition<Type,MeshType>> boundaryConditions_;

protected:

    // Allocate data
    void allocate();

    // Transfer data from other mesh field
    void transferData(meshField<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("meshField");

    // Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    // List type
    typedef PtrList<meshLevel<Type,MeshType>> listType;

    // Main constructor: construct from name and mesh

    meshField
    (
        const word& name,
        const fvMesh& fvMsh,
        const IOobject::readOption r = NO_READ,
        const IOobject::writeOption w = NO_WRITE,
        const bool registerObject = false,
        const bool initBCs = false
    );

    // Copy constructors: construct from field or tmp of field with option to
    // clone or transfer the boundary conditions

    meshField
    (
        const meshField<Type,MeshType>& field,
        const bool registerObject = false,
        const bool copyBCs = false
    );

    meshField
    (
        const word& name,
        const meshField<Type,MeshType>& field,
        const bool registerObject = false,
        const bool copyBCs = false
    );

    meshField
    (
        const tmp<meshField<Type,MeshType>>& field,
        const bool registerObject = false,
        const bool copyBCs = false
    );

    meshField
    (
        const word& name,
        const tmp<meshField<Type,MeshType>>& field,
        const bool registerObject = false,
        const bool copyBCs = false
    );

    // Clone
    inline autoPtr<meshField<Type,MeshType>> clone() const
    {
        return autoPtr<meshField<Type,MeshType>>
        (
            new meshField<Type,MeshType>(*this)
        );
    }

    // Destructor
    ~meshField();

    // Access to the levels. Use the operator() for this because it doesn't
    // perform any bounds checking which is expensive.
    inline meshLevel<Type,MeshType>& operator[](label l)
    {
        return *const_cast<meshLevel<Type,MeshType>*>
        (
            listType::operator()(l)
        );
    }

    // Const access to the levels. Use the operator() for this because it
    // doesn't perform any bounds checking which is expensive.
    inline const meshLevel<Type,MeshType>& operator[](label l) const
    {
        return *listType::operator()(l);
    }

    // Number of levels
    inline label size() const { return listType::size(); }

    // Initialize ghost cells to zero
    void initGhosts();

    // Initialize ghost cells to the given value
    void initGhosts(const Type&);

    // Return the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Old time field pointer
    inline meshField<Type,MeshType>* oldTimePtr() const
    {
        return oldTimePtr_;
    }

    // Access to the old time field. Return the current time of no old time
    // field is set
    inline meshField<Type,MeshType>& oldTime()
    {
        return oldTimePtr_ == nullptr ? *this : *oldTimePtr_;
    }

    inline const meshField<Type,MeshType>& oldTime() const
    {
        return oldTimePtr_ == nullptr ? *this : *oldTimePtr_;
    }

    // Access to the boundary conditions
    inline PtrList<boundaryCondition<Type,MeshType>>&
    boundaryConditions()
    {
        return boundaryConditions_;
    }

    // Const access to the boundary conditions
    inline const PtrList<boundaryCondition<Type,MeshType>>&
    boundaryConditions() const
    {
        return boundaryConditions_;
    }

    // Add boundary conditions
    void addBoundaryConditions();

    // Correct boundary conditions
    void correctBoundaryConditions(const bool = false);

    // Check if boundary conditions are singular, which is true if all boundary
    // conditions are empty, parallel, periodic or Neumann
    bool singularBoundaryConditions();

    // Set the old time field
    void setOldTime();

    // Write data (not supported)
    bool writeData(Ostream& os) const { return os; }

    // Member operators

    void operator=(const meshField<Type,MeshType>&);
    void operator=(const tmp<meshField<Type,MeshType>>&);
    void operator=(const Type&);
    void operator=(const zero);

    void operator+=(const meshField<Type,MeshType>&);
    void operator+=(const tmp<meshField<Type,MeshType>>&);

    void operator-=(const meshField<Type,MeshType>&);
    void operator-=(const tmp<meshField<Type,MeshType>>&);

    void operator*=(const meshField<scalar,MeshType>&);
    void operator*=(const tmp<meshField<scalar,MeshType>>&);

    void operator/=(const meshField<scalar,MeshType>&);
    void operator/=(const tmp<meshField<scalar,MeshType>>&);

    void operator+=(const Type&);
    void operator-=(const Type&);

    void operator*=(const scalar&);
    void operator/=(const scalar&);

    // Templated functions are needed for types which can be cast to another
    // type

    template<class Type2> void operator=(const meshField<Type2,MeshType>&);
    template<class Type2> void operator=(const tmp<meshField<Type2,MeshType>>&);

    template<class Type2> void operator+=(const meshField<Type2,MeshType>&);
    template<class Type2> void operator+=(const tmp<meshField<Type2,MeshType>>&);

    template<class Type2> void operator-=(const meshField<Type2,MeshType>&);
    template<class Type2> void operator-=(const tmp<meshField<Type2,MeshType>>&);

    template<class Type2> void operator=(const Type2&);
    template<class Type2> void operator+=(const Type2&);
    template<class Type2> void operator-=(const Type2&);

    // Allow a parent linear system to make changes to private data
    template<class SType,class,class>
    friend class linearSystem;
};

}

}

}

#include "meshFieldFunctions.H"
#include "meshFieldStencilFunctions.H"

#ifdef NoRepository
    #include "meshField.C"
#endif

#endif
