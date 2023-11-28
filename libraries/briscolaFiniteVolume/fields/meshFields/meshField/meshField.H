#ifndef meshField_H
#define meshField_H

#include "block.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "PtrList.H"

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
class restrictionScheme;

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

    // Boundary conditions, assigned in the order:
    //
    //  1) non-empty boundaries,
    //  2) parallel/periodic faces,
    //  3) parallel/periodic edges,
    //  4) parallel/periodic vertices and
    //  5) empty boundaries.
    //
    //  This order is important, so that parallel boundary conditions can copy
    //  normal boundaries and empty boundaries can copy both normal and
    //  parallel/periodic boundaries.
    PtrList<boundaryCondition<Type,MeshType>> boundaryConditions_;

    // Restriction scheme associated with this field
    autoPtr<restrictionScheme<Type,MeshType>> reScheme_;

protected:

    // Allocate data and specify depth
    void allocate(const bool);

    // Set field pointers in levels, needed when data is transferred
    void setFieldPointers();

    // Transfer data from other mesh field
    void transfer(meshField<Type,MeshType>&);

    // Check if the given field's depth is compatible with this field's depth to
    // execute a member operator
    template<class Type2>
    inline void checkMemberOperatorArgDepth
    (
        const meshField<Type2,MeshType>& F
    ) const
    {
        if (this->deep() && F.shallow())
            FatalErrorInFunction
                << "Member operator with a shallow field argument on a deep "
                << "field is not supported." << endl << abort(FatalError);
    }

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
        const bool initBCs = false,
        const bool deep = false
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
        #ifdef FULLDEBUG

        if (l < 0 || l >= listType::size())
            FatalErrorInFunction
                << "Index " << l << " out of range 0 ... " << listType::size()-1
                << endl << abort(FatalError);

        #endif

        return *const_cast<meshLevel<Type,MeshType>*>
        (
            listType::operator()(l)
        );
    }

    // Const access to the levels. Use the operator() for this because it
    // doesn't perform any bounds checking which is expensive.
    inline const meshLevel<Type,MeshType>& operator[](label l) const
    {
        #ifdef FULLDEBUG

        if (l < 0 || l >= listType::size())
            FatalErrorInFunction
                << "Index " << l << " out of range 0 ... " << listType::size()-1
                << endl << abort(FatalError);

        #endif

        return *listType::operator()(l);
    }

    // Alias access function
    inline meshLevel<Type,MeshType>& level(const int l)
    {
        return this->operator[](l);
    }

    // Alias access function
    inline const meshLevel<Type,MeshType>& level(const int l) const
    {
        return this->operator[](l);
    }

    // Access to the first level
    inline meshLevel<Type,MeshType>& mshPart()
    {
        return this->operator[](0);
    }

    // Const access to the first level
    inline const meshLevel<Type,MeshType>& mshPart() const
    {
        return this->operator[](0);
    }

    // Access to the first direction of the first level
    inline meshDirection<Type,MeshType>& direction()
    {
        return this->operator[](0)[0];
    }

    // Const access to the first direction of the first level
    inline const meshDirection<Type,MeshType>& direction() const
    {
        return this->operator[](0)[0];
    }

    // Access to the given direction of the first level
    inline meshDirection<Type,MeshType>& direction(const int d)
    {
        return this->operator[](0)[d];
    }

    // Const access to the given direction of the first level
    inline const meshDirection<Type,MeshType>& direction(const int d) const
    {
        return this->operator[](0)[d];
    }

    // Access to the given direction and level
    inline meshDirection<Type,MeshType>& direction
    (
        const int l,
        const int d
    )
    {
        return this->operator[](l)[d];
    }

    // Const access to the given direction and level
    inline const meshDirection<Type,MeshType>& direction
    (
        const int l,
        const int d
    ) const
    {
        return this->operator[](l)[d];
    }

    // Number of levels
    inline label size() const { return listType::size(); }

    // Alias to size
    inline label numberOfLevels() const { return size(); }

    // Number of directions
    inline label numberOfDirections() const
    {
        return MeshType::numberOfDirections;
    }

    // Return the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Old time field pointer
    inline meshField<Type,MeshType>* oldTimePtr() const
    {
        return oldTimePtr_;
    }

    // Access to the old time field. Return the current time if no old time
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
    void correctBoundaryConditions();

    // Correct parallel boundary conditions
    void correctParallelBoundaryConditions();

    // Correct periodic boundary conditions
    void correctPeriodicBoundaryConditions();

    // Correct communicating (parallel and periodic) boundary conditions
    void correctCommBoundaryConditions();

    // Correct non-communicating boundary conditions
    void correctNonCommBoundaryConditions();

    // Set the old time field
    void setOldTime();

    // Make the field deep
    void makeDeep();

    // Make the field shallow
    void makeShallow();

    // Make the field deep or shallow
    inline void make(const bool deep)
    {
        if (deep)
            makeDeep();
        else
            makeShallow();
    }

    // The field is deep if it contains more than one level
    inline bool deep() const
    {
        return listType::size() > 1;
    }

    // Shallow if not deep
    inline bool shallow() const
    {
        return !deep();
    }

    // Set the restriction operator that is associated with this field
    void setRestrictionScheme(const word scheme);

    // Restrict the field to lower levels. If no restriction scheme is set, set
    // the default restriction scheme first. A shallow field will be made deep.
    void restrict();

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
    template<class Type2> void operator+=
    (
        const tmp<meshField<Type2,MeshType>>&
    );

    template<class Type2> void operator-=(const meshField<Type2,MeshType>&);
    template<class Type2> void operator-=
    (
        const tmp<meshField<Type2,MeshType>>&
    );

    template<class Type2> void operator=(const Type2&);
    template<class Type2> void operator+=(const Type2&);
    template<class Type2> void operator-=(const Type2&);

    // Allow a parent linear system to make changes to private data
    template<class SType,class,class>
    friend class linearSystem;

    #include "meshLevelFirstDirectionFunctions.H"
    #include "meshLevelDirectionFunctions.H"
    #include "meshFieldDirectionFunctions.H"
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
