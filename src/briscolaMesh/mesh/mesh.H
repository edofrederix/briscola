#ifndef mesh_H
#define mesh_H

#include "decomposition.H"
#include "geometry.H"
#include "part.H"
#include "faceLabel.H"
#include "edgeLabel.H"
#include "vertexLabel.H"

// Forward declaration of MPI communicator
struct ompi_communicator_t;
typedef ompi_communicator_t* MPI_Comm;

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

    // Mask to keep track of set boundaries
    labelBlock setBoundaryMask_;

    // Bool to check if mesh is structured or unstructured. If the brick
    // topology is structured, so is the mesh.
    bool structured_;

    // Check if mesh directions are rectilinear
    labelVector rectilinear_;

    // Check if mesh directions are uniform
    labelVector uniform_;

    // Bounding box
    faceScalar boundingBox_;

    // Distributed communication graph
    MPI_Comm comm_;

    // Add a boundary to the list of boundaries
    void addBoundary
    (
        const word name,
        const word type,
        const labelTensor T,
        const labelVector offset,
        const vector mappedOffset = vector::zero,
        const labelVector neighborOffset = zeroXYZ,
        const label neighborProcNum = -1
    );

    // Generate brick-internal boundaries
    void generateBrickInternalBoundaries();

    // Generate brick-external boundary
    void generateBrickExternalBoundaries();

    // Generate link boundaries from brick link
    void generateLinkBoundaries(const brickLink&);

    // Generate patch boundaries
    void generatePatchBoundaries();

    // Generate boundaries
    void generateBoundaries();

    // Set boundary mask
    void setBoundaryMask();

    // Reorder boundaries in the order:
    //
    //  1) patch/mapped boundaries,
    //  2) parallel/periodic faces,
    //  3) parallel/periodic edges,
    //  4) parallel/periodic vertices and
    //  5) empty boundaries.
    //
    //  This order is important, so that parallel boundary conditions can copy
    //  patch/mapped boundaries and empty boundaries can copy both patch/mapped
    //  and parallel/periodic boundaries.
    void reorderBoundaries();

    // Set communication tags in parallel and periodic boundaries. Tags avoid
    // inconsistent communication.
    void setCommTags();

    // Set boundary extension labels
    void setBoundaryExtension();

    // Set the list of empty boundary offsets
    void setEmptyBoundaryOffsets();

    // Set the distributed communication graph
    void setDistributedCommGraph();

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

    // Get boundary from offset vector
    inline const boundary& b(const labelVector& bo) const
    {
        forAll(boundaries_, i)
            if (boundaries_[i].offset() == bo)
                return boundaries_[i];

        FatalErrorInFunction
            << "Invalid boundary offset" << endl << abort(FatalError);

        return boundaries_[0];
    }

    // Return the set boundary mask
    inline const labelBlock& setBoundaryMask() const
    {
        return setBoundaryMask_;
    }

    // Structured mesh bool
    inline bool structured() const { return structured_; }

    // Return rectilinear bools
    inline labelVector rectilinear() const { return rectilinear_; }

    // Return uniform bools
    inline labelVector uniform() const { return uniform_; }

    // Bounding box
    inline faceScalar boundingBox() const { return boundingBox_; }

    // Access to the distributed communicator graph
    const MPI_Comm& comm() const { return comm_; }

    // Find the cell index of the given point at given mesh level. Return
    // -unitXYZ if it is not found. Tree search, but should be overloaded by
    // mesh specializations for better performance.
    virtual labelVector findCell(const vector& p, const label l = 0) const;

    // Find the cell indices of the given list of points at given mesh level
    virtual List<labelVector> findCells(const vectorList&, const label) const;

    // Given a part size, compute its coarsened child size subject to this
    // mesh's properties
    labelVector coarsen(const labelVector) const;

    // Disallow default bitwise assignment
    void operator=(const mesh&);
};

}

}

#endif
