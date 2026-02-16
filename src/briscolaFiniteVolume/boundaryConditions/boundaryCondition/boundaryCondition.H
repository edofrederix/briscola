#ifndef boundaryCondition_H
#define boundaryCondition_H

#include "boundaryConditionBaseTypes.H"

#include "block.H"

#include "faceVector.H"
#include "faceScalar.H"
#include "stencil.H"
#include "diagStencil.H"

#include "boundary.H"

#include "PtrList.H"
#include "NamedEnum.H"
#include "products.H"

#include "meshFieldFwd.H"
#include "faceFieldFwd.H"
#include "linearSystemFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh;

template<class Type, class MeshType>
class boundaryCondition
:
    public refCount
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the mesh level
    meshLevel<Type,MeshType>& level_;

    // Level number to which this boundary condition belongs
    label l_;

    // Reference to the base boundary
    const boundary& b_;

    // Copy of the dictionary
    dictionary dict_;

public:

    // Runtime type information
    TypeName("boundaryCondition");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        boundaryCondition,
        dictionary,
        (
            const meshLevel<Type,MeshType>& level,
            const boundary& b
        ),
        (level, b)
    );

    // Construct from components
    boundaryCondition
    (
        const meshLevel<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    boundaryCondition(const boundaryCondition<Type,MeshType>&);

    // Copy constructor for an alternative level
    boundaryCondition
    (
        const boundaryCondition<Type,MeshType>&,
        const meshLevel<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const = 0;

    // Construct and return a clone for alternative level
    virtual tmp<boundaryCondition<Type,MeshType>> clone
    (
        const meshLevel<Type,MeshType>&
    ) const = 0;

    // Selector which determines type based on the provided boundary type
    static autoPtr<boundaryCondition<Type,MeshType>> New
    (
        const meshLevel<Type,MeshType>&,
        const boundary&
    );

    // Check if cast to derived type is possible
    template<class BoundaryType>
    inline bool castable() const
    {
        return dynamic_cast<const BoundaryType*>(this);
    }

    // Cast to derived type, if possible. Returns error otherwise.
    template<class BoundaryType>
    inline BoundaryType& cast()
    {
        BoundaryType* bPtr = dynamic_cast<BoundaryType*>(this);

        if (!bPtr)
        {
            FatalErrorInFunction
                << "Boundary condition cannot be cast to "
                << BoundaryType::typeName
                << " because it is " << this->type()
                << abort(FatalError);
        }

        return *bPtr;
    }

    // Cast to const derived type reference
    template<class BoundaryType>
    inline const BoundaryType& cast() const
    {
        return
            const_cast<boundaryCondition<Type,MeshType>*>(this)
          ->cast<BoundaryType>();
    }

    // Determine a level's boundary condition base type for a global boundary
    // offset. Works only on rectilinear meshes. Gives error for inconsistent
    // boundaries.
    static boundaryConditionBaseType globalBaseType
    (
        const meshLevel<Type,MeshType>&,
        const labelVector
    );

    // Destructor
    virtual ~boundaryCondition();

    // Return the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the mesh level
    inline meshField<Type,MeshType>& level()
    {
        return level_;
    }

    // Const access to the mesh level
    inline const meshField<Type,MeshType>& level() const
    {
        return level_;
    }

    // Level number
    inline label levelNum() const { return l_; }

    // Access to the base boundary
    inline const boundary& b() const { return b_; }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Return the boundary offset
    inline labelVector offset() const
    {
        return b_.offset();
    }

    // Boundary offset degree
    inline label offsetDegree() const
    {
        return b_.offsetDegree();
    }

    // Return the start point of the internal cell boundary
    labelVector S(const label d) const;

    // Return the end point of the internal cell boundary
    labelVector E(const label d) const;

    // Return the shape of the internal cell boundary
    labelVector N(const label d) const;

    // Access to the face centers of the field this BC belongs to
    const faceField<vector,MeshType>& faceCenters() const;

    // Access to the face normals of the field this BC belongs to
    const faceField<vector,MeshType>& faceNormals() const;

    // Access to the face areas of the field this BC belongs to
    const faceField<scalar,MeshType>& faceAreas() const;

    // Access to the cell volumes of the field this BC belongs to
    const meshField<scalar,MeshType>& cellVolumes() const;

    // Access to the cell centers of the field this BC belongs to
    const meshField<vector,MeshType>& cellCenters() const;

    // Access to the face deltas of the field this BC belongs to
    const faceField<scalar,MeshType>& faceDeltas() const;

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const = 0;

    // Eliminated boundary condition
    virtual bool eliminated() const = 0;

    // Eliminate ghosts from the linear system at this BC's level. Cannot be
    // templated because of virtual functions.

    virtual void eliminateGhosts
    (
        linearSystem<stencil,Type,MeshType>&
    );

    virtual void eliminateGhosts
    (
        linearSystem<diagStencil,Type,MeshType>&
    );

    // Eliminate ghosts from the linear system at this BC's level and direction.
    // Cannot be templated because of virtual functions.

    virtual void eliminateGhosts
    (
        linearSystem<stencil,Type,MeshType>&,
        const label
    );

    virtual void eliminateGhosts
    (
        linearSystem<diagStencil,Type,MeshType>&,
        const label
    );

    // Prepare the evaluation of the boundary condition at this BC's level and
    // optionally at given direction
    virtual void prepare();
    virtual void prepare(const label);

    // Evaluate the boundary condition at this BC's level and optionally at
    // given direction
    virtual void evaluate();
    virtual void evaluate(const label);
};

// Helper function for global boundary condition type

template<class Type, class MeshType>
boundaryConditionBaseType globalBoundaryConditionBaseType
(
    const meshLevel<Type,MeshType>& level,
    const labelVector bo
)
{
    return boundaryCondition<Type,MeshType>::globalBaseType(level,bo);
}

}

}

}

#ifdef NoRepository
    #include "boundaryCondition.C"
#endif

#endif
