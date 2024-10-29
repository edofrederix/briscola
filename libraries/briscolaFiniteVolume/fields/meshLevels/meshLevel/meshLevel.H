#ifndef meshLevel_H
#define meshLevel_H

#include "PtrList.H"
#include "List.H"

#include "block.H"
#include "Vector.H"

#include "meshDirection.H"

// #define BOUNDARYEXCHANGE 1

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh;

template<class Type, class MeshType>
class meshField;

template<class Type, class MeshType>
class meshLevel
:
    public PtrList<meshDirection<Type,MeshType>>,
    public refCount
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Level number
    label l_;

    // Pointer to the mesh field, if any
    meshField<Type,MeshType>* mshFieldPtr_;

    // Allocate data
    void allocate();

    // Set level pointers in directions, needed when data is transferred
    void setLevelPointers();

    // Transfer data from other mesh level
    void transfer(meshLevel<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("meshLevel");

    // Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    // List type
    typedef PtrList<meshDirection<Type,MeshType>> listType;

    // Main constructor: construct from owner field, mesh and level number

    meshLevel
    (
        meshField<Type,MeshType>&,
        const fvMesh&,
        const label
    );

    // Copy constructors: construct from mesh level, tmp or mesh

    meshLevel(const meshLevel<Type,MeshType>&);
    meshLevel(const meshLevel<Type,MeshType>&, const zero&);
    meshLevel(const meshLevel<Type,MeshType>&, const Type&);
    meshLevel(const meshLevel<Type,MeshType>&, const List<Type>&);

    meshLevel(const tmp<meshLevel<Type,MeshType>>&);
    meshLevel(const tmp<meshLevel<Type,MeshType>>&, const zero&);
    meshLevel(const tmp<meshLevel<Type,MeshType>>&, const Type&);
    meshLevel(const tmp<meshLevel<Type,MeshType>>&, const List<Type>&);

    meshLevel(const fvMesh&, const label);
    meshLevel(const fvMesh&, const label, const zero&);
    meshLevel(const fvMesh&, const label, const Type&);
    meshLevel(const fvMesh&, const label, const List<Type>&);

    // Clone

    inline autoPtr<meshLevel<Type,MeshType>> clone() const
    {
        return autoPtr<meshLevel<Type,MeshType>>
        (
            new meshLevel<Type,MeshType>(*this)
        );
    }

    inline autoPtr<meshLevel<Type,MeshType>> clone(const zero&) const
    {
        return autoPtr<meshLevel<Type,MeshType>>
        (
            new meshLevel<Type,MeshType>(*this, Zero)
        );
    }

    inline autoPtr<meshLevel<Type,MeshType>> clone(const Type& v) const
    {
        return autoPtr<meshLevel<Type,MeshType>>
        (
            new meshLevel<Type,MeshType>(*this, v)
        );
    }

    inline autoPtr<meshLevel<Type,MeshType>> clone(const List<Type>& v) const
    {
        return autoPtr<meshLevel<Type,MeshType>>
        (
            new meshLevel<Type,MeshType>(*this, v)
        );
    }

    // Destructor
    virtual ~meshLevel();

    // Access to the directions. Use the operator() for this because it doesn't
    // perform any bounds checking which is expensive.
    inline meshDirection<Type,MeshType>& operator[](label d)
    {
        #ifdef FULLDEBUG

        if (d < 0 || d >= listType::size())
            FatalErrorInFunction
                << "Index " << d << " out of range 0 ... " << listType::size()-1
                << endl << abort(FatalError);

        #endif

        return *const_cast<meshDirection<Type,MeshType>*>
        (
            listType::operator()(d)
        );
    }

    // Const access to the directions. Use the operator() for this because it
    // doesn't perform any bounds checking which is expensive.
    inline const meshDirection<Type,MeshType>& operator[](label d) const
    {
        #ifdef FULLDEBUG

        if (d < 0 || d >= listType::size())
            FatalErrorInFunction
                << "Index " << d << " out of range 0 ... " << listType::size()-1
                << endl << abort(FatalError);

        #endif

        return *listType::operator()(d);
    }

    // Alias access function
    inline meshDirection<Type,MeshType>& direction(const int d)
    {
        return this->operator[](d);
    }

    // Alias const access function
    inline const meshDirection<Type,MeshType>& direction(const int d) const
    {
        return this->operator[](d);
    }

    // Access the first direction
    inline meshDirection<Type,MeshType>& direction()
    {
        return this->operator[](0);
    }

    // Const access the first direction
    inline const meshDirection<Type,MeshType>& direction() const
    {
        return this->operator[](0);
    }

    // Number of directions
    inline label size() const { return MeshType::numberOfDirections; }

    // Alias to size
    inline label numberOfDirections() const { return size(); }

    // Return the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Pointer to the mesh field
    inline const meshField<Type,MeshType>* mshFieldPtr() const
    {
        return mshFieldPtr_;
    }

    // Const reference to the mesh field (unchecked!)
    inline const meshField<Type,MeshType>& mshField() const
    {
        return *mshFieldPtr_;
    }

    // Access to the mesh field (unchecked!)
    inline meshField<Type,MeshType>& mshField()
    {
        return *mshFieldPtr_;
    }

    // Level number
    inline label levelNum() const { return l_; }

    // Mesh part
    inline const part& mshPart() const { return fvMsh_[l_]; }

    // Add boundary conditions only of not yet added
    inline void addBoundaryConditions()
    {
        if (mshFieldPtr_ && mshFieldPtr_->boundaryConditions().size() == 0)
        {
            mshFieldPtr_->addBoundaryConditions();
        }
    }

    // Add immersed boundary conditions only of not yet added
    inline void addImmersedBoundaryConditions()
    {
        if
        (
            mshFieldPtr_
         && mshFieldPtr_->immersedBoundaryConditions().size() == 0
        )
        {
            mshFieldPtr_->addImmersedBoundaryConditions();
        }
    }

    // Correct boundary conditions
    void correctBoundaryConditions();

    // Allow boundary conditions to be corrected on a const level
    void correctBoundaryConditions() const
    {
        const_cast<meshLevel<Type,MeshType>*>(this)
            ->correctBoundaryConditions();
    }

    // Correct unset boundary conditions as homogeneous Neumann
    void correctUnsetBoundaryConditions();

    // Correct domain boundary conditions
    void correctDomainBoundaryConditions();

    // Correct empty boundary conditions
    void correctEmptyBoundaryConditions();

    // Correct communicating boundary conditions
    void correctCommsBoundaryConditions();

    // Correct parallel boundary conditions
    void correctParallelBoundaryConditions();

    // Correct periodic boundary conditions
    void correctPeriodicBoundaryConditions();

    // Correct non-eliminated boundary conditions
    void correctNonEliminatedBoundaryConditions();

    // Correct eliminated boundary conditions
    void correctEliminatedBoundaryConditions();

    // Correct immersed boundary conditions
    void correctImmersedBoundaryConditions();

    // Member operators

    void operator=(const meshLevel<Type,MeshType>&);
    void operator=(const tmp<meshLevel<Type,MeshType>>&);
    void operator=(const Type&);
    void operator=(const List<Type>&);
    void operator=(const zero);

    void operator+=(const meshLevel<Type,MeshType>&);
    void operator+=(const tmp<meshLevel<Type,MeshType>>&);

    void operator-=(const meshLevel<Type,MeshType>&);
    void operator-=(const tmp<meshLevel<Type,MeshType>>&);

    void operator*=(const meshLevel<scalar,MeshType>&);
    void operator*=(const tmp<meshLevel<scalar,MeshType>>&);

    void operator/=(const meshLevel<scalar,MeshType>&);
    void operator/=(const tmp<meshLevel<scalar,MeshType>>&);

    void operator+=(const Type&);
    void operator+=(const List<Type>&);

    void operator-=(const Type&);
    void operator-=(const List<Type>&);

    void operator*=(const scalar&);
    void operator*=(const List<scalar>&);

    void operator/=(const scalar&);
    void operator/=(const List<scalar>&);

    // Templated functions are needed for types which can be cast to another
    // type

    #define TYPE2LEVEL meshLevel<Type2,MeshType>

    template<class Type2> void operator=(const TYPE2LEVEL&);
    template<class Type2> void operator=(const tmp<TYPE2LEVEL>&);

    template<class Type2> void operator+=(const TYPE2LEVEL&);
    template<class Type2> void operator+=(const tmp<TYPE2LEVEL>&);

    template<class Type2> void operator-=(const TYPE2LEVEL&);
    template<class Type2> void operator-=(const tmp<TYPE2LEVEL>&);

    template<class Type2> void operator*=(const TYPE2LEVEL&);
    template<class Type2> void operator*=(const tmp<TYPE2LEVEL>&);

    template<class Type2> void operator/=(const TYPE2LEVEL&);
    template<class Type2> void operator/=(const tmp<TYPE2LEVEL>&);

    #undef TYPE2LEVEL

    template<class Type2> void operator=(const Type2&);
    template<class Type2> void operator=(const List<Type2>&);

    template<class Type2> void operator+=(const Type2&);
    template<class Type2> void operator+=(const List<Type2>&);

    template<class Type2> void operator-=(const Type2&);
    template<class Type2> void operator-=(const List<Type2>&);

    template<class Type2> void operator*=(const Type2&);
    template<class Type2> void operator*=(const List<Type2>&);

    template<class Type2> void operator/=(const Type2&);
    template<class Type2> void operator/=(const List<Type2>&);

    // Allow the parent field to make changes to private data
    friend class meshField<Type,MeshType>;

    #include "meshLevelFirstDirectionFunctions.H"
    #include "meshLevelDirectionFunctions.H"
};

}

}

}

#include "meshLevelFunctions.H"
#include "meshLevelStencilFunctions.H"

#ifdef NoRepository
    #include "meshLevel.C"
#endif

#endif
