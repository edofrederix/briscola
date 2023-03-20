#ifndef boundaryCondition_H
#define boundaryCondition_H

#include "boundaryConditionBaseTypes.H"

#include "faceVector.H"
#include "faceScalar.H"

#include "partPatch.H"

#include "PtrList.H"
#include "NamedEnum.H"

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
class boundaryCondition
:
    public refCount
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the mesh field
    meshField<Type,MeshType>& mshField_;

    // Reference to the part patch
    const partPatch& patch_;

    // Copy of the dictionary
    dictionary dict_;

    // Disallow default bitwise assignment
    void operator=(const boundaryCondition<Type,MeshType>&);

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
            const meshField<Type,MeshType>& mshField,
            const partPatch& patch
        ),
        (mshField, patch)
    );

    // Construct from components
    boundaryCondition
    (
        const meshField<Type,MeshType>&,
        const partPatch&
    );

    // Copy constructor
    boundaryCondition(const boundaryCondition<Type,MeshType>&);

    // Copy constructor resetting field reference
    boundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const = 0;

    // Construct and return a clone setting field reference
    virtual tmp<boundaryCondition<Type,MeshType>> clone
    (
        const meshField<Type,MeshType>& field
    ) const = 0;

    // Selector for given type
    static autoPtr<boundaryCondition<Type,MeshType>> New
    (
        const meshField<Type,MeshType>&,
        const partPatch&,
        const word
    );

    // Selector for a normal boundary condition
    static autoPtr<boundaryCondition<Type,MeshType>> NewBoundary
    (
        const meshField<Type,MeshType>&,
        const partPatch&
    );

    // Selector for a periodic boundary condition
    static autoPtr<boundaryCondition<Type,MeshType>> NewPeriodic
    (
        const meshField<Type,MeshType>&,
        const partPatch& p
    );

    // Selector for a parallel boundary condition
    static autoPtr<boundaryCondition<Type,MeshType>> NewParallel
    (
        const meshField<Type,MeshType>&,
        const partPatch& p
    );

    // Determine a field's boundary condition base type for a global boundary
    // offset. Works only on rectilinear meshes. Gives error for inconistent
    // boundaries.
    static boundaryConditionBaseType globalBaseType
    (
        const meshField<Type,MeshType>&,
        const labelVector
    );

    // Destructor
    virtual ~boundaryCondition();

    // Return the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the mesh field
    inline meshField<Type,MeshType>& mshField()
    {
        return mshField_;
    }

    // Const access to the mesh field
    inline const meshField<Type,MeshType>& mshField() const
    {
        return mshField_;
    }

    // Access to the part patch
    inline const partPatch& patch() const { return patch_; }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Return the boundary offset
    inline labelVector boundaryOffset() const
    {
        return patch_.boundaryOffset();
    }

    // Boundary offset degree
    inline label boundaryOffsetDegree() const
    {
        return patch_.boundaryOffsetDegree();
    }

    // Return the transformation tensor
    inline labelTensor T() const
    {
        return patch_.T();
    }

    // Return master bool for shifted boundaries
    virtual bool master() const
    {
        return patch_.master();
    }

    // Return slave bool for shifted boundaries
    virtual bool slave() const
    {
        return patch_.slave();
    }

    // Access to the face centers of the field this BC belongs to
    const meshField<faceVector,MeshType>& faceCenters() const;

    // Access to the face normals of the field this BC belongs to
    const meshField<faceVector,MeshType>& faceNormals() const;

    // Access to the face areas of the field this BC belongs to
    const meshField<faceScalar,MeshType>& faceAreas() const;

    // Access to the cell volumes of the field this BC belongs to
    const meshField<scalar,MeshType>& cellVolumes() const;

    // Access to the cell centers of the field this BC belongs to
    const meshField<vector,MeshType>& cellCenters() const;

    // Access to the face deltas of the field this BC belongs to
    const meshField<faceScalar,MeshType>& faceDeltas() const;

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const = 0;

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label) = 0;

    // Evaluate the boundary condition for given level. Optionally force
    // homogeneous BC
    virtual void evaluate(const label, const bool = false) = 0;
};

// Helper function for global boundary condition type

template<class Type, class MeshType>
boundaryConditionBaseType globalBoundaryConditionBaseType
(
    const meshField<Type,MeshType>& field,
    const labelVector bo
)
{
    return boundaryCondition<Type,MeshType>::globalBaseType(field,bo);
}

}

}

}

#ifdef NoRepository
    #include "boundaryCondition.C"
#endif

#endif
