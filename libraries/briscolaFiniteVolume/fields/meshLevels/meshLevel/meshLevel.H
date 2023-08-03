#ifndef meshLevel_H
#define meshLevel_H

#include "PtrList.H"
#include "List.H"

#include "block.H"
#include "Vector.H"

#include "meshDirection.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh;

// template<class Type, class MeshType>
// class boundaryCondition;

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

    // Transfer data from other mesh level
    void transferData(meshLevel<Type,MeshType>&);

    // Update active cells based on boundary conditions
    void updateActiveCells();

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

    // Destructor
    virtual ~meshLevel();

    // Access to the directions. Use the operator() for this because it doesn't
    // perform any bounds checking which is expensive.
    inline meshDirection<Type,MeshType>& operator[](label d)
    {
        return *const_cast<meshDirection<Type,MeshType>*>
        (
            listType::operator()(d)
        );
    }

    // Const access to the directions. Use the operator() for this because it
    // doesn't perform any bounds checking which is expensive.
    inline const meshDirection<Type,MeshType>& operator[](label d) const
    {
        return *listType::operator()(d);
    }

    // Number of directions
    inline label size() const { return listType::size(); }

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

    // Part level
    inline const partLevel& level() const { return fvMsh_[l_]; }

    // Add boundary conditions only of not yet added
    inline void addBoundaryConditions()
    {
        if (mshFieldPtr_ && mshFieldPtr_->boundaryConditions().size() == 0)
        {
            mshFieldPtr_->addBoundaryConditions();
        }
    }

    // Correct boundary conditions for this level
    void correctBoundaryConditions(const bool = false);

    // Correct boundary boundary conditions for this level
    void correctBoundaryBoundaryConditions();

    // Correct parallel boundary conditions for this level
    void correctParallelBoundaryConditions();

    // Correct periodic boundary conditions for this level
    void correctPeriodicBoundaryConditions();

    // Correct communicating (parallel and periodic) boundary conditions for
    // this level
    void correctCommBoundaryConditions();

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

    template<class Type2> void operator=(const meshLevel<Type2,MeshType>&);
    template<class Type2> void operator=(const tmp<meshLevel<Type2,MeshType>>&);

    template<class Type2> void operator+=(const meshLevel<Type2,MeshType>&);
    template<class Type2> void operator+=(const tmp<meshLevel<Type2,MeshType>>&);

    template<class Type2> void operator-=(const meshLevel<Type2,MeshType>&);
    template<class Type2> void operator-=(const tmp<meshLevel<Type2,MeshType>>&);

    template<class Type2> void operator=(const Type2&);
    template<class Type2> void operator=(const List<Type2>&);
    template<class Type2> void operator+=(const Type2&);
    template<class Type2> void operator+=(const List<Type2>&);
    template<class Type2> void operator-=(const Type2&);
    template<class Type2> void operator-=(const List<Type2>&);

    // Allow the parent field to make changes to private data
    friend class meshField<Type,MeshType>;
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
