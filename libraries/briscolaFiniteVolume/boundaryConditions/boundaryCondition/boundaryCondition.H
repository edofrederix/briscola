#ifndef boundaryCondition_H
#define boundaryCondition_H

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

// Enumeration listing the base boundary condition types
enum boundaryConditionBaseType
{
    DUMMYBC,
    EMPTYBC,
    PARALLELBC,
    PERIODICBC,
    DIRICHLETBC,
    NEUMANNBC,
    ROBINBC
};

const NamedEnum<boundaryConditionBaseType,7> boundaryConditionBaseTypeNames;

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
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new boundaryCondition<Type,MeshType>(*this)
        );
    }

    // Construct and return a clone setting field reference
    virtual tmp<boundaryCondition<Type,MeshType>> clone
    (
        const meshField<Type,MeshType>& field
    ) const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new boundaryCondition<Type,MeshType>(field,*this)
        );
    }

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

    // Return master bool
    inline bool master() const
    {
        return patch_.master();
    }

    // Return slave bool
    inline bool slave() const
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
    virtual boundaryConditionBaseType baseType() const
    {
        NotImplemented;
        return DUMMYBC;
    }

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label)
    {
        NotImplemented;
    }

    // Evaluate the boundary condition for given level. Optionally force
    // homogeneous BC
    virtual void evaluate(const label, const bool = false)
    {
        NotImplemented;
    }

    // Boundary conditions on shifted boundaries may be non-eliminated. For
    // example, for Dirichlet-like boundary conditions the value on the boundary
    // is known. Also for slave parallel/periodic boundary conditions the value
    // is known. By default, boundary conditions are eliminated.
    virtual bool eliminated() const
    {
        return true;
    }
};

}

}

}

#define makeBoundaryCondition(Type,MeshType)                                    \
                                                                                \
    typedef boundaryCondition<Type,MeshType>                                    \
        boundaryCondition##Type##MeshType;                                      \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        boundaryCondition##Type##MeshType,                                      \
        "boundaryCondition",                                                    \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        boundaryCondition##Type##MeshType,                                      \
        dictionary                                                              \
    );                                                                          \
                                                                                \
    template class boundaryCondition<Type,MeshType>;

#define makeBoundaryConditionType(BC, Type, MeshType)                           \
                                                                                \
    typedef BC##BoundaryCondition<Type,MeshType>                                \
        BC##BoundaryCondition##Type##MeshType;                                  \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        BC##BoundaryCondition##Type##MeshType,                                  \
        #BC,                                                                    \
        0                                                                       \
    );                                                                          \
                                                                                \
    boundaryCondition<Type,MeshType>::                                          \
        adddictionaryConstructorToTable<BC##BoundaryCondition##Type##MeshType>  \
            add##BC##BoundaryCondition##Type##MeshType##ConstructorToTable_;    \
                                                                                \
    template class BC##BoundaryCondition<Type,MeshType>;

#endif
