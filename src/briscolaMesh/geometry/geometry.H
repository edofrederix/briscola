#ifndef geometry_H
#define geometry_H

#include "meshData.H"
#include "brick.H"
#include "brickTopology.H"
#include "patch.H"
#include "patchPair.H"

namespace Foam
{

namespace briscola
{

class geometry
:
    public meshData
{
    // List of bricks
    PtrList<brick> bricks_;

    // Brick connectivity
    autoPtr<brickTopology> topology_;

    // List of patches
    PtrList<patch> patches_;

    // List of patch pairs
    PtrList<patchPair> patchPairs_;

    // Create bricks
    void createBricks();

    // Create the brick connectivity
    void createBrickTopology();

    // Create patches
    void createPatches();

    // Create the default patch (requires topology to be set)
    void createDefaultPatch();

    /// Add given patch and possibly split up the patch if it has faces that
    /// share the same brick
    void addPatch(const patch& p);

    // Create patch pairs
    void createPatchPairs();

    // Check consistency of the geometry
    void checkPatchConsistency() const;

    // Align bricks w.r.t. the first in the topology map. Only possible for
    // structured topologies.
    void alignBricks();

public:

    // Construct from IOdictionary
    geometry(const IOdictionary&);

    // Copy constructor. Copies bricks and re-generates everything else.
    geometry(const geometry& geo);

    // Destructor
    ~geometry();

    // Reference to the list of bricks
    PtrList<brick>& bricks()
    {
        return bricks_;
    }

    // Const reference to the list of bricks
    const PtrList<brick>& bricks() const
    {
        return bricks_;
    }

    // Const reference to the brick topology
    const brickTopology& topology() const
    {
        return topology_();
    }

    // List of patches
    const PtrList<patch>& patches() const
    {
        return patches_;
    }

    // List of patch pairs
    const PtrList<patchPair>& patchPairs() const
    {
        return patchPairs_;
    }
};

}

}

#endif
