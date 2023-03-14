#ifndef mesh_H
#define mesh_H

#include "decomposition.H"
#include "geometry.H"
#include "partLevel.H"

namespace Foam
{

namespace briscola
{

class mesh
:
    public geometry,
    public PtrList<partLevel>
{
    // Mesh decomposition
    autoPtr<decomposition> decomp_;

    // This processor's part patches
    PtrList<partPatch> partPatches_;

    // Lower and upper master patch vectors. Indicate for each patch if it is a
    // master patch.
    labelVector lowerMasterPatch_;
    labelVector upperMasterPatch_;

    // Lower and upper boundary patch vectors. Indicate for each patch if it is
    // a boundary patch.
    labelVector lowerBoundaryPatch_;
    labelVector upperBoundaryPatch_;

    // Bool to check if mesh is structured or unstructured. If the brick
    // topology is structured, so is the mesh.
    bool structured_;

    // Check if mesh directions are rectilinear
    labelVector rectilinear_;

    // Check if mesh directions are uniform
    labelVector uniform_;

    // Add a part patch to the list of part patches
    void addPartPatch
    (
        const word name,
        const word type,
        const labelTensor T,
        const labelVector boundaryOffset,
        const labelVector neighborOffset = zeroXYZ,
        const label neighborProcNum = -1
    );

    // Generate brick-internal part patches
    void generateBrickInternalPartPatches();

    // Generate brick-external part patches
    void generateBrickExternalPartPatches();

    // Generate link patches from brick link
    void generateLinkPartPatches(const brickLink&);

    // Generate boundary patches
    void generateBoundaryPartPatches();

    // Generate part patches
    void generatePartPatches();

    // Set communication tags in parallel and periodic patches. Tags avoid
    // inconsistent communication.
    void setCommTags();

    // Generate part levels
    void generatePartLevels();

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

    // Access to the part patches
    inline const PtrList<partPatch>& partPatches() const
    {
        return partPatches_;
    }

    // Access to the lower master patch vector
    inline labelVector lowerMasterPatch() const
    {
        return lowerMasterPatch_;
    }

    // Access to the lower slave patch vector
    inline labelVector lowerSlavePatch() const
    {
        return unitXYZ - lowerMasterPatch_;
    }

    // Access to the upper master patch vector
    inline labelVector upperMasterPatch() const
    {
        return upperMasterPatch_;
    }

    // Access to the upper slave patch vector
    inline labelVector upperSlavePatch() const
    {
        return unitXYZ - upperMasterPatch_;
    }

    // Access to the lower boundary patch vector
    inline labelVector lowerBoundaryPatch() const
    {
        return lowerBoundaryPatch_;
    }

    // Access to the upper boundary patch vector
    inline labelVector upperBoundaryPatch() const
    {
        return upperBoundaryPatch_;
    }

    // Structured mesh bool
    inline bool structured() const { return structured_; }

    // Return rectilinear bools
    inline labelVector rectilinear() const { return rectilinear_; }

    // Return uniform bools
    inline labelVector uniform() const { return uniform_; }

    //- Disallow default bitwise assignment
    void operator=(const mesh&);
};

}

}

#endif
