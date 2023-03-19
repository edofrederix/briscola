#ifndef meshDirection_H
#define meshDirection_H

#include "PtrList.H"
#include "List.H"

#include "block.H"
#include "Vector.H"

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
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Level number
    label l_;

    // Direction number
    label d_;

    // Internal cell block start
    labelVector S_;

    // Internal cell block end
    labelVector E_;

    // Internal cell block size
    labelVector N_;

    // Active cell block start
    labelVector Sa_;

    // Active cell block end
    labelVector Ea_;

    // Active cell block size
    labelVector Na_;

    // Pointer to the mesh level, if any
    meshLevel<Type,MeshType>* mshLevelPtr_;

    // Set internal cells
    void setInternalCells();

    // Update active cells based on boundary conditions
    void updateActiveCells();

    // Allocate data
    void allocate(const labelVector);

    // Transfer data from other mesh direction
    void transferData(meshDirection<Type,MeshType>&);

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
        meshLevel<Type,MeshType>&,
        const fvMesh&,
        const label,
        const label
    );

    // Copy constructors: construct from mesh direction, tmp or mesh

    meshDirection(const meshDirection<Type,MeshType>&);
    meshDirection(const meshDirection<Type,MeshType>&, const zero&);
    meshDirection(const meshDirection<Type,MeshType>&, const Type&);

    meshDirection(const tmp<meshDirection<Type,MeshType>>&);
    meshDirection(const tmp<meshDirection<Type,MeshType>>&, const zero&);
    meshDirection(const tmp<meshDirection<Type,MeshType>>&, const Type&);

    meshDirection(const fvMesh&, const label, const label);
    meshDirection(const fvMesh&, const label, const label, const zero&);
    meshDirection(const fvMesh&, const label, const label, const Type&);

    // Destructor
    virtual ~meshDirection();

    // Return the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Mesh level pointer
    inline const meshLevel<Type,MeshType>* mshLevelPtr() const
    {
        return mshLevelPtr_;
    }

    // Mesh level (unchecked!)
    inline const meshLevel<Type,MeshType>& mshLevel() const
    {
        return *mshLevelPtr_;
    }

    // Level number
    inline label levelNum() const { return l_; }

    // Direction number
    inline label directionNum() const { return d_; }

    // Part level
    inline const partLevel& level() const { return fvMsh_[l_]; }

    // Access to base block
    inline block<Type>& B() { return *this; }

    // Const access to the base
    inline const block<Type>& B() const { return *this; }

    // Return the internal cell block start
    inline labelVector S() const { return S_; }

    // Return the internal cell block end
    inline labelVector E() const { return E_; }

    // Return the internal cell block size
    inline labelVector N() const { return N_; }

    // Return the active cell block start
    inline labelVector Sa() const { return Sa_; }

    // Return the active cell block end
    inline labelVector Ea() const { return Ea_; }

    // Return the active cell block size
    inline labelVector Na() const { return Na_; }

    // Internal cell boundaries

    // Return the start point of the internal cell boundary at given offset
    inline labelVector boundaryStart(const labelVector& bo) const
    {
        return labelVector
        (
            bo.x() == 0 ? 0 : bo.x() == -1 ? S_.x() : E_.x()-1,
            bo.y() == 0 ? 0 : bo.y() == -1 ? S_.y() : E_.y()-1,
            bo.z() == 0 ? 0 : bo.z() == -1 ? S_.z() : E_.z()-1
        );
    }

    // Return the shape of the internal cell boundary at given offset
    inline labelVector boundaryN(const labelVector& bo) const
    {
        return labelVector
        (
            bo.x() == 0 ? N_.x() : 1,
            bo.y() == 0 ? N_.y() : 1,
            bo.z() == 0 ? N_.z() : 1
        );
    }

    // Return the end of the internal cell boundary at given offset
    inline labelVector boundaryEnd(const labelVector& bo) const
    {
        return boundaryStart(bo) + boundaryN(bo);
    }

    // Return the start point of the active cell boundary at given offset
    inline labelVector activeBoundaryStart(const labelVector& bo) const
    {
        return labelVector
        (
            bo.x() == 0 ? 0 : bo.x() == -1 ? Sa_.x() : Ea_.x()-1,
            bo.y() == 0 ? 0 : bo.y() == -1 ? Sa_.y() : Ea_.y()-1,
            bo.z() == 0 ? 0 : bo.z() == -1 ? Sa_.z() : Ea_.z()-1
        );
    }

    // Return the shape of the active cell boundary at given offset
    inline labelVector activeBoundaryN(const labelVector& bo) const
    {
        return labelVector
        (
            bo.x() == 0 ? Na_.x() : 1,
            bo.y() == 0 ? Na_.y() : 1,
            bo.z() == 0 ? Na_.z() : 1
        );
    }

    // Return the end of the active cell boundary at given offset
    inline labelVector activeBoundaryEnd(const labelVector& bo) const
    {
        return activeBoundaryStart(bo) + activeBoundaryN(bo);
    }

    // Return true if the boundary at given offset is shifted
    inline bool shifted(const labelVector& bo) const
    {
        for (int dim = 0; dim < 3; dim++)
            if (bo[dim] != 0 && MeshType::shift[d_][dim] != 0.0)
                return true;

        return false;
    }

    // Return the shape without ghost cells as a diagonal tensor
    inline labelTensor shapeTensor() const
    {
        return labelTensor(N_.x(),0,0,0,N_.y(),0,0,0,N_.z());
    }

    // Return the number of internal cells
    inline int size() const { return cmptProduct(N_); }

    // Return the number of active cells
    inline int activeSize() const { return cmptProduct(Na_); }

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

    // (i,j,k) element access operator. Account for ghost cell layers
    inline Type& operator()
    (
        const int i,
        const int j,
        const int k
    )
    {
        return blockType::operator()(i+1,j+1,k+1);
    }

    // (i,j,k) const element access operator. Account for ghost cell layers
    inline const Type& operator()
    (
        const int i,
        const int j,
        const int k
    ) const
    {
        return blockType::operator()(i+1,j+1,k+1);
    }

    // (i,j,k) element access operator. Account for ghost cell layers
    inline Type& operator() (const labelVector& ijk)
    {
        return blockType::operator()(ijk+unitXYZ);
    }

    // (i,j,k) const element access operator. Account for ghost cell layers
    inline const Type& operator() (const labelVector& ijk) const
    {
        return blockType::operator()(ijk+unitXYZ);
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

    template<class Type2> void operator=(const meshDirection<Type2,MeshType>&);
    template<class Type2> void operator=(const tmp<meshDirection<Type2,MeshType>>&);

    template<class Type2> void operator+=(const meshDirection<Type2,MeshType>&);
    template<class Type2> void operator+=(const tmp<meshDirection<Type2,MeshType>>&);

    template<class Type2> void operator-=(const meshDirection<Type2,MeshType>&);
    template<class Type2> void operator-=(const tmp<meshDirection<Type2,MeshType>>&);

    template<class Type2> void operator=(const Type2&);
    template<class Type2> void operator+=(const Type2&);
    template<class Type2> void operator-=(const Type2&);

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
