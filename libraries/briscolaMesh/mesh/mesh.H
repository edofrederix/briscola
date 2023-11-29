#ifndef mesh_H
#define mesh_H

#include "decomposition.H"
#include "geometry.H"
#include "part.H"
#include "faceLabel.H"
#include "edgeLabel.H"
#include "vertexLabel.H"

namespace Foam
{

namespace briscola
{

class mesh
:
    public geometry,
    public PtrList<part>
{
    // Mesh decomposition
    autoPtr<decomposition> decomp_;

    // This processor's boundaries
    PtrList<boundary> boundaries_;

    // Boundary master face labels, indicating for each processor which face
    // boundaries are master. Orientation is processor-specific.
    List<faceLabel> faceBoundaryMasterPerProc_;
    List<edgeLabel> edgeBoundaryMasterPerProc_;
    List<vertexLabel> vertexBoundaryMasterPerProc_;

    // Boundary type face labels, indicating for each processor the type of the
    // face boundaries. Orientation is processor-specific.
    List<faceLabel> faceBoundaryTypePerProc_;
    List<edgeLabel> edgeBoundaryTypePerProc_;
    List<vertexLabel> vertexBoundaryTypePerProc_;

    // Patch extension labels
    edgeLabel edgePatchExtension_;
    vertexLabel vertexPatchExtension_;

    // List of empty boundary offsets, which are boundary offsets in which no
    // boundary is present
    List<labelVector> emptyPatchOffsets_;

    // Bool to check if mesh is structured or unstructured. If the brick
    // topology is structured, so is the mesh.
    bool structured_;

    // Check if mesh directions are rectilinear
    labelVector rectilinear_;

    // Check if mesh directions are uniform
    labelVector uniform_;

    // Bounding box
    faceScalar boundingBox_;

    // Add a boundary to the list of boundaries
    void addBoundary
    (
        const word name,
        const word type,
        const labelTensor T,
        const labelVector offset,
        const labelVector neighborOffset = zeroXYZ,
        const label neighborProcNum = -1
    );

    // Generate brick-internal boundaries
    void generateBrickInternalBoundaries();

    // Generate brick-external boundary
    void generateBrickExternalBoundaries();

    // Generate link boundaries from brick link
    void generateLinkBoundaries(const brickLink&);

    // Generate domain boundaries
    void generateDomainBoundaries();

    // Generate boundaries
    void generateBoundaries();

    // Set communication tags in parallel and periodic boundaries. Tags avoid
    // inconsistent communication.
    void setCommTags();

    // Set boundary master and type label lists
    void setPatchLabels();

    // Set boundary extension labels
    void setPatchExtension();

    // Generate the list of empty boundary offsets
    void setEmptyPatchOffsets();

    // Generate mesh part levels
    void generateLevels();

public:

    // Runtime type information
    TypeName("mesh");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        mesh,
        dictionary,
        (autoPtr<mesh>& mshPtr),
        (mshPtr)
    );

    // Construct from IOdictionary
    explicit mesh(const IOdictionary&);

    // Copy constructor
    mesh(const mesh& msh);

    // Transfer constructor and reset arg
    mesh(autoPtr<mesh>&);

    // Transfer or copy constructor
    mesh(mesh& msh, bool reuse);

    // Clone mesh
    virtual autoPtr<mesh> clone() const
    {
        return autoPtr<mesh>(new mesh(*this));
    }

    // Selector
    static autoPtr<mesh> New(const IOdictionary&);

    // Destructor
    virtual ~mesh();

    // Check if cast to derived type is possible
    template<class MeshType>
    inline bool castable() const
    {
        return dynamic_cast<const MeshType*>(this);
    }

    // Cast to derived type, if possible. Returns error otherwise.
    template<class MeshType>
    inline const MeshType& cast() const
    {
        const MeshType* mshPtr = dynamic_cast<const MeshType*>(this);

        if (!mshPtr)
        {
            FatalErrorInFunction
                << "Mesh cannot be cast to " << MeshType::typeName
                << " because it is " << this->type()
                << abort(FatalError);
        }

        return *mshPtr;
    }

    // Mesh decomposition
    inline const decomposition& decomp() const
    {
        return decomp_;
    }

    // Access to the boundaries
    inline const PtrList<boundary>& boundaries() const
    {
        return boundaries_;
    }

    // Const access to the boundary master label lists

    inline const List<faceLabel>& faceBoundaryMasterPerProc() const
    {
        return faceBoundaryMasterPerProc_;
    }

    inline const List<edgeLabel>& edgeBoundaryMasterPerProc() const
    {
        return edgeBoundaryMasterPerProc_;
    }

    inline const List<vertexLabel>& vertexBoundaryMasterPerProc() const
    {
        return vertexBoundaryMasterPerProc_;
    }

    // Const access to the boundary type label list

    inline const List<faceLabel>& faceBoundaryTypePerProc() const
    {
        return faceBoundaryTypePerProc_;
    }

    inline const List<edgeLabel>& edgeBoundaryTypePerProc() const
    {
        return edgeBoundaryTypePerProc_;
    }

    inline const List<vertexLabel>& vertexBoundaryTypePerProc() const
    {
        return vertexBoundaryTypePerProc_;
    }

    // Return this processor's boundary master labels

    inline const faceLabel& faceBoundaryMaster() const
    {
        return faceBoundaryMasterPerProc_[Pstream::myProcNo()];
    }

    inline const edgeLabel& edgeBoundaryMaster() const
    {
        return edgeBoundaryMasterPerProc_[Pstream::myProcNo()];
    }

    inline const vertexLabel& vertexBoundaryMaster() const
    {
        return vertexBoundaryMasterPerProc_[Pstream::myProcNo()];
    }

    // Return this processor's boundary slave face label

    inline faceLabel facePatchSlave() const
    {
        return
            pTraits<faceLabel>::one
          - faceBoundaryMasterPerProc_[Pstream::myProcNo()];
    }

    inline edgeLabel edgePatchSlave() const
    {
        return
            pTraits<edgeLabel>::one
          - edgeBoundaryMasterPerProc_[Pstream::myProcNo()];
    }

    inline vertexLabel vertexPatchSlave() const
    {
        return
            pTraits<vertexLabel>::one
          - vertexBoundaryMasterPerProc_[Pstream::myProcNo()];
    }

    // Return this processor's boundary type face label

    inline const faceLabel& faceBoundaryType() const
    {
        return faceBoundaryTypePerProc_[Pstream::myProcNo()];
    }

    inline const edgeLabel& edgeBoundaryType() const
    {
        return edgeBoundaryTypePerProc_[Pstream::myProcNo()];
    }

    inline const vertexLabel& vertexBoundaryType() const
    {
        return vertexBoundaryTypePerProc_[Pstream::myProcNo()];
    }

    // Return boundary extension

    inline const edgeLabel& edgePatchExtension() const
    {
        return edgePatchExtension_;
    }

    inline const vertexLabel& vertexPatchExtension() const
    {
        return vertexPatchExtension_;
    }

    // Return the list of empty boundary offsets
    inline const List<labelVector>& emptyPatchOffsets() const
    {
        return emptyPatchOffsets_;
    }

    // Structured mesh bool
    inline bool structured() const { return structured_; }

    // Return rectilinear bools
    inline labelVector rectilinear() const { return rectilinear_; }

    // Return uniform bools
    inline labelVector uniform() const { return uniform_; }

    // Bounding box
    inline faceScalar boundingBox() const { return boundingBox_; }

    // Find the cell index of the given point at given mesh level. Return
    // -unitXYZ if it is not found. Tree search, but should be overloaded by
    // mesh specializations for better performance.
    virtual labelVector findCell(const vector& p, const label l = 0) const;

    // Find the cell indices of the given list of points at given mesh level
    virtual List<labelVector> findCells(const vectorList&, const label) const;

    // Disallow default bitwise assignment
    void operator=(const mesh&);
};

}

}

#endif
