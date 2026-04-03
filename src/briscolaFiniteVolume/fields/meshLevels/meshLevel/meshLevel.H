#ifndef meshLevel_H
#define meshLevel_H

#include "meshFieldFwd.H"
#include "FastPtrList.H"
#include "List.H"

#include "block.H"
#include "Vector.H"

#include "boundaryCondition.H"
#include "meshDirection.H"

#define BOUNDARYEXCHANGE 1

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh;

template<class Type, class MeshType>
class boundaryExchange;

template<class Type, class MeshType>
class meshLevel
:
    public FastPtrList<meshDirection<Type,MeshType>>,
    public refCount
{
public:

    // Number of ghost cell layers
    static const label ghosts = meshDirection<Type,MeshType>::ghosts;

private:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Level number
    label l_;

    // Pointer to the mesh field, if any
    meshField<Type,MeshType>* fieldPtr_;

    // Boundary conditions
    FastPtrList<boundaryCondition<Type,MeshType>> boundaryConditions_;

    // Boundary exchange (defined on field level to enable collective operations
    // on all levels)
    autoPtr<boundaryExchange<Type,MeshType>> bExchangePtr_;

    // Singular boundary condition bool (set when adding boundary conditions)
    bool singular_;

    // Allocate data
    void allocate();

    // Set level pointers in directions, needed when data is transferred
    void setLevelPointers();

    // Transfer data from other mesh level
    void transfer(meshLevel<Type,MeshType>&);

    // Check if the given mesh level matches our level
    template<class Type2>
    void checkLevel(const meshLevel<Type2,MeshType>& L) const
    {
        if (L.levelNum() != l_)
            FatalErrorInFunction
                << "Invalid level" << endl << abort(FatalError);
    }

public:

    // Runtime type information
    TypeName("meshLevel");

    // Data type
    typedef Type dataType;

    // Mesh type
    typedef MeshType meshType;

    // Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    // List type
    typedef FastPtrList<meshDirection<Type,MeshType>> listType;

    // Main constructor: construct from owner field, mesh and level number
    meshLevel(meshField<Type,MeshType>& field, const label l);

    // Copy constructors from level
    meshLevel(const meshLevel<Type,MeshType>& L);
    meshLevel(const tmp<meshLevel<Type,MeshType>>&);

    // Standalone constructors: construct without parent field
    meshLevel(const fvMesh& fvMsh, const label l);
    meshLevel(const fvMesh& fvMsh, const label l, const zero&);
    meshLevel(const fvMesh& fvMsh, const label l, const Type& v);
    meshLevel(const fvMesh& fvMsh, const label l, const List<Type>& v);

    // Temporary level constructors

    static tmp<meshLevel<Type,MeshType>> New
    (
        meshField<Type,MeshType>& field,
        const label l
    )
    {
        return
            tmp<meshLevel<Type,MeshType>>
            (
                new meshLevel<Type,MeshType>(field,l)
            );
    }

    static tmp<meshLevel<Type,MeshType>> New
    (
        const fvMesh& fvMsh,
        const label l
    )
    {
        return
            tmp<meshLevel<Type,MeshType>>
            (
                new meshLevel<Type,MeshType>(fvMsh,l)
            );
    }

    // Clone
    inline autoPtr<meshLevel<Type,MeshType>> clone() const
    {
        return autoPtr<meshLevel<Type,MeshType>>
        (
            new meshLevel<Type,MeshType>(*this)
        );
    }

    // Destructor
    virtual ~meshLevel();

    // Access to the directions
    inline meshDirection<Type,MeshType>& operator[](label d)
    {
        return listType::operator[](d);
    }

    // Const access to the directions
    inline const meshDirection<Type,MeshType>& operator[](label d) const
    {
        return listType::operator[](d);
    }

    // Alias access function
    inline meshDirection<Type,MeshType>& direction(const label d)
    {
        return listType::operator[](d);
    }

    // Alias const access function
    inline const meshDirection<Type,MeshType>& direction(const label d) const
    {
        return listType::operator[](d);
    }

    // Access the first direction
    inline meshDirection<Type,MeshType>& direction()
    {
        return listType::operator[](0);
    }

    // Const access the first direction
    inline const meshDirection<Type,MeshType>& direction() const
    {
        return listType::operator[](0);
    }

    // Number of directions
    inline label size() const { return MeshType::numberOfDirections; }

    // Alias to size
    inline label numberOfDirections() const { return size(); }

    // Return the finite volume mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return the mesh
    inline const mesh& msh() const { return fvMsh_.msh(); }

    // Pointer to the mesh field
    inline const meshField<Type,MeshType>* fieldPtr() const
    {
        return fieldPtr_;
    }

    // Const reference to the mesh field (unchecked!)
    inline const meshField<Type,MeshType>& field() const
    {
        return *fieldPtr_;
    }

    // Access to the mesh field (unchecked!)
    inline meshField<Type,MeshType>& field()
    {
        return *fieldPtr_;
    }

    // Level number
    inline label levelNum() const { return l_; }

    // Return the mesh level
    inline const level& lvl() const { return fvMsh_[l_]; }

    // Access to the boundary conditions
    inline FastPtrList<boundaryCondition<Type,MeshType>>& boundaryConditions()
    {
        return boundaryConditions_;
    }

    // Const access to the boundary conditions
    inline const FastPtrList<boundaryCondition<Type,MeshType>>&
    boundaryConditions() const
    {
        return boundaryConditions_;
    }

    // Get boundary condition from offset vector
    inline boundaryCondition<Type,MeshType>& boundary
    (
        const labelVector& bo
    )
    {
        forAll(boundaryConditions_, i)
            if (boundaryConditions_[i].offset() == bo)
                return boundaryConditions_[i];

        FatalErrorInFunction
            << "Invalid boundary condition offset" << endl << abort(FatalError);

        return boundaryConditions_[0];
    }

    // Get const reference to boundary condition from offset vector
    inline const boundaryCondition<Type,MeshType>& boundary
    (
        const labelVector& bo
    ) const
    {
        return const_cast<meshField<Type,MeshType>*>(this)->boundary(bo);
    }

    // Add boundary conditions only of not yet added
    void addBoundaryConditions();

    // Add immersed boundary conditions only of not yet added
    inline void addImmersedBoundaryConditions()
    {
        if
        (
            fieldPtr_
         && fieldPtr_->immersedBoundaryConditions().size() == 0
        )
        {
            fieldPtr_->addImmersedBoundaryConditions();
        }
    }

    // Transfer boundary conditions from another level
    void transferBoundaryConditions(const meshLevel<Type,MeshType>&);

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
    void correctPatchBoundaryConditions();

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

    // Correct agglomeration data
    void correctAggData();

    // Check if the boundary conditions make this level singular
    bool singular() const
    {
        const_cast<meshLevel<Type,MeshType>*>(this)->addBoundaryConditions();

        return singular_;
    }

    // Return a component
    tmp<meshLevel<cmptType,MeshType>> component(const label) const;

    // Set a component from a directional list of values
    void replace(const label, const List<cmptType>&);

    // Set a component from a component level
    void replace(const label, const meshLevel<cmptType,MeshType>&);
    void replace(const label, const tmp<meshLevel<cmptType,MeshType>>&);

    // Set maximum/minimum of level and given value
    void max(const Type&);
    void min(const Type&);

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

    // Befriend parallel exchange
    template<class, class>
    friend class boundaryExchange;

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
