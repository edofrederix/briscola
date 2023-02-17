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

    // Mesh sizes of this processor's brick's part
    labelVector N_;

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

    // Declare name of the class and its debug switch
    ClassName("mesh");

    // Construct from IOdictionary
    explicit mesh(const IOdictionary&);

    //- Disallow default bitwise copy construction
    mesh(const mesh&);

    // Destructor
    ~mesh();

    // Mesh decomposition
    inline const decomposition& decomp() const
    {
        return decomp_;
    }

    // Mesh sizes of this processor's brick's part
    inline const labelVector& N() const
    {
        return N_;
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

    //- Disallow default bitwise assignment
    void operator=(const mesh&);
};

}

}

#endif
