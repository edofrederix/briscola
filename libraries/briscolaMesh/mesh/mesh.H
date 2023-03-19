#ifndef mesh_H
#define mesh_H

#include "decomposition.H"
#include "geometry.H"
#include "partLevel.H"
#include "faceLabel.H"

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

    // Patch master face labels, indicating for each processor which face part
    // patches are master. Orientation is processor-specific.
    List<faceLabel> patchMasterPerProc_;

    // Patch type face labels, indicating for each processor the type of the
    // face part patches. Orientation is processor-specific.
    List<faceLabel> patchTypePerProc_;

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

    // Const access to the patch master face label list
    inline const List<faceLabel> patchMasterPerProc() const
    {
        return patchMasterPerProc_;
    }

    // Const access to the patch type face label list
    inline const List<faceLabel> patchTypePerProc() const
    {
        return patchTypePerProc_;
    }

    // Return this processor's patch master face label
    inline faceLabel patchMaster() const
    {
        return patchMasterPerProc_[Pstream::myProcNo()];
    }

    // Return this processor's patch slave face label
    inline faceLabel patchSlave() const
    {
        return
            pTraits<faceLabel>::one
          - patchMasterPerProc_[Pstream::myProcNo()];
    }

    // Return this processor's patch type face label
    inline faceLabel patchType() const
    {
        return patchTypePerProc_[Pstream::myProcNo()];
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
