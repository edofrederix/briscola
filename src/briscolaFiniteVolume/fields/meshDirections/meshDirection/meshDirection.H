#ifndef meshDirection_H
#define meshDirection_H

#include "PtrList.H"
#include "List.H"

#include "block.H"
#include "Vector.H"
#include "faceLabel.H"

#include "fvMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class meshField;

template<class Type, class MeshType>
class meshLevel;

template<class Type, class MeshType>
class meshDirection
:
    // Keep the block private, to avoid unwarranted modification or inconsistent
    // indexing. Direct access to the base block remains possible with the B()
    // member function.

    private block<Type>
{
public:

    // Number of ghost cell layers
    static const label ghosts = 1;

private:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Level number
    label l_;

    // Direction number
    label d_;

    // Store internal cell block start and end indices for faster access
    faceLabel I_;

    // Pointer to the mesh level, if any
    meshLevel<Type,MeshType>* mshLevelPtr_;

    // Allocate data
    void allocate();

    // Transfer data from other mesh direction
    void transfer(meshDirection<Type,MeshType>& D);

public:

    // Runtime type information
    TypeName("meshDirection");

    // Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    // Block type
    typedef block<Type> blockType;

    // Provide public access to the block's refCount member functions
    using blockType::refCount::count;
    using blockType::refCount::unique;
    using blockType::refCount::operator--;
    using blockType::refCount::operator++;

    // Main constructor: construct from owner level, mesh, level number and
    // direction number

    meshDirection
    (
        meshLevel<Type,MeshType>& mshLevel,
        const fvMesh& fvMsh,
        const label l,
        const label d
    );

    // Copy constructors: construct from mesh direction, tmp or mesh

    meshDirection
    (
        const meshDirection<Type,MeshType>& D
    );

    meshDirection
    (
        const meshDirection<Type,MeshType>& D,
        const zero&
    );

    meshDirection
    (
        const meshDirection<Type,MeshType>& D,
        const Type&
    );

    meshDirection
    (
        const meshDirection<Type,MeshType>& D,
        const List<Type>&
    );

    meshDirection
    (
        const tmp<meshDirection<Type,MeshType>>& tD
    );

    meshDirection
    (
        const tmp<meshDirection<Type,MeshType>>& tD,
        const zero&
    );

    meshDirection
    (
        const tmp<meshDirection<Type,MeshType>>& tD,
        const Type&
    );

    meshDirection
    (
        const tmp<meshDirection<Type,MeshType>>& tD,
        const List<Type>&
    );

    // Standalone constructors: construct without parent level with optional
    // ghost layers

    meshDirection
    (
        const fvMesh& fvMsh,
        const label l,
        const label d
    );

    meshDirection
    (
        const fvMesh& fvMsh,
        const label l,
        const label d,
        const zero& z
    );

    meshDirection
    (
        const fvMesh& fvMsh,
        const label l,
        const label d,
        const Type& v
    );

    meshDirection
    (
        const fvMesh& fvMsh,
        const label l,
        const label d,
        const List<Type>& v
    );

    // Clone

    inline autoPtr<meshDirection<Type,MeshType>> clone() const
    {
        return autoPtr<meshDirection<Type,MeshType>>
        (
            new meshDirection<Type,MeshType>(*this)
        );
    }

    inline autoPtr<meshDirection<Type,MeshType>> clone(const zero&) const
    {
        return autoPtr<meshDirection<Type,MeshType>>
        (
            new meshDirection<Type,MeshType>(*this, Zero)
        );
    }

    inline autoPtr<meshDirection<Type,MeshType>> clone(const Type& v) const
    {
        return autoPtr<meshDirection<Type,MeshType>>
        (
            new meshDirection<Type,MeshType>(*this, v)
        );
    }

    inline autoPtr<meshDirection<Type,MeshType>> clone
    (
        const List<Type>& v
    ) const
    {
        return autoPtr<meshDirection<Type,MeshType>>
        (
            new meshDirection<Type,MeshType>(*this, v)
        );
    }

    // Destructor
    virtual ~meshDirection();

    // Return the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Mesh level pointer
    inline const meshLevel<Type,MeshType>* mshLevelPtr() const
    {
        return mshLevelPtr_;
    }

    // Const reference to the mesh level (unchecked!)
    inline const meshLevel<Type,MeshType>& mshLevel() const
    {
        return *mshLevelPtr_;
    }

    // Access to the mesh level (unchecked!)
    inline meshLevel<Type,MeshType>& mshLevel()
    {
        return *mshLevelPtr_;
    }

    // Level number
    inline label levelNum() const { return l_; }

    // Direction number
    inline label directionNum() const { return d_; }

    // Mesh part
    inline const part& mshPart() const { return fvMsh_[l_]; }

    // Access to base block
    inline block<Type>& B() { return *this; }

    // Const access to the base
    inline const block<Type>& B() const { return *this; }

    // Internal cell start and end indices
    inline const faceLabel& I() const { return I_; }

    // Internal cell block size
    inline labelVector N() const { return I_.upper()-I_.lower(); }

    // Return the number of internal cells
    inline int size() const { return cmptProduct(N()); }

    // Linear element access operator will not work efficiently due to
    // ghost cells, so let's prevent its use
    inline Type& operator()(const int i)
    {
        NotImplemented;
        return this->B()(i);
    }

    inline const Type& operator()(const int i) const
    {
        NotImplemented;
        return this->B()(i);
    }

    // (i,j,k) element access operator. Account for ghost cell layers.
    inline Type& operator()
    (
        const int i,
        const int j,
        const int k
    )
    {
        return blockType::operator()(i+ghosts, j+ghosts, k+ghosts);
    }

    // (i,j,k) const element access operator. Account for ghost cell layers.
    inline const Type& operator()
    (
        const int i,
        const int j,
        const int k
    ) const
    {
        return blockType::operator()(i+ghosts, j+ghosts, k+ghosts);
    }

    // (i,j,k) element access operator. Account for ghost cell layers.
    inline Type& operator() (const labelVector& ijk)
    {
        return blockType::operator()(ijk + ghosts*unitXYZ);
    }

    // (i,j,k) const element access operator. Account for ghost cell layers.
    inline const Type& operator() (const labelVector& ijk) const
    {
        return blockType::operator()(ijk + ghosts*unitXYZ);
    }

    // Member operators

    void operator=(const meshDirection<Type,MeshType>&);
    void operator=(const tmp<meshDirection<Type,MeshType>>&);
    void operator=(const Type&);
    void operator=(const zero);

    void operator+=(const meshDirection<Type,MeshType>&);
    void operator+=(const tmp<meshDirection<Type,MeshType>>&);

    void operator-=(const meshDirection<Type,MeshType>&);
    void operator-=(const tmp<meshDirection<Type,MeshType>>&);

    void operator*=(const meshDirection<scalar,MeshType>&);
    void operator*=(const tmp<meshDirection<scalar,MeshType>>&);

    void operator/=(const meshDirection<scalar,MeshType>&);
    void operator/=(const tmp<meshDirection<scalar,MeshType>>&);

    void operator+=(const Type&);
    void operator-=(const Type&);

    void operator*=(const scalar&);
    void operator/=(const scalar&);

    // Templated functions are needed for types which can be cast to another
    // type

    #define TYPE2DIR meshDirection<Type2,MeshType>

    template<class Type2> void operator=(const TYPE2DIR&);
    template<class Type2> void operator=(const tmp<TYPE2DIR>&);

    template<class Type2> void operator+=(const TYPE2DIR&);
    template<class Type2> void operator+=(const tmp<TYPE2DIR>&);

    template<class Type2> void operator-=(const TYPE2DIR&);
    template<class Type2> void operator-=(const tmp<TYPE2DIR>&);

    template<class Type2> void operator*=(const TYPE2DIR&);
    template<class Type2> void operator*=(const tmp<TYPE2DIR>&);

    template<class Type2> void operator/=(const TYPE2DIR&);
    template<class Type2> void operator/=(const tmp<TYPE2DIR>&);

    #undef TYPE2DIR

    template<class Type2> void operator=(const Type2&);
    template<class Type2> void operator+=(const Type2&);
    template<class Type2> void operator-=(const Type2&);
    template<class Type2> void operator*=(const Type2&);
    template<class Type2> void operator/=(const Type2&);

    // Allow the parent level to make changes to private data
    friend class meshLevel<Type,MeshType>;
};

}

}

}

#include "meshDirectionFunctions.H"
#include "meshDirectionStencilFunctions.H"

#ifdef NoRepository
    #include "meshDirection.C"
#endif

#endif
