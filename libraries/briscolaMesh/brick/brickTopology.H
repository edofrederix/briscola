#ifndef brickTopology_H
#define brickTopology_H

#include "brick.H"
#include "PtrList.H"
#include "brickLink.H"
#include "brickLinks.H"
#include "brickTopologyMap.H"
#include "labelList.H"
#include "labelVector.H"

namespace Foam
{

namespace briscola
{

class geometry;
class patch;

class brickTopology
{
    // Reference to the geometry object
    const geometry& geo_;

    // For each brick a collection of its links to other bricks
    PtrList<brickLinks> links_;

    // Bool to check if bricks are structured or unstructured
    bool structured_;

    // Bool to check if the bricks are aligned
    bool aligned_;

    // Bool to check if the bricks are on a rectilinear grid
    bool rectilinear_;

    // Brick topology map, only set if the bricks are structured. Its local
    // coordinate system is aligned with that of the first brick.
    brickTopologyMap map_;

    // Set face links
    void setFaceLinks();

    // Set periodic face links
    void setPeriodicFaceLinks();

    // Set edge and vertex links through face link paths
    void setEdgeAndVertexLinks();

    // Recursive function to build the brick topology map. Returns false when
    // failed, meaning that bricks are unstructured
    bool buildTopologyMap
    (
        const brick&,
        const labelTensor,
        labelVector&,
        labelList&,
        labelBlock&
    ) const;

public:

    // Construct from geometry
    brickTopology(const geometry&);

    // Destructor
    virtual ~brickTopology();

    // Find the shortest path between two brick numbers via faces. Option to
    // ignore periodic faces.
    labelList shortestFacePath
    (
        const label&,
        const label&,
        const bool = false
    ) const;

    // Return brick links
    inline const PtrList<brickLinks>& links() const
    {
        return links_;
    }

    // Structured brick topology bool
    inline bool structured() const { return structured_; }

    // Aligned brick topology bool
    inline bool aligned() const { return aligned_; }

    // Rectilinear brick topology bool
    inline bool rectilinear() const { return rectilinear_; }

    // Reference to the brick topology map, which is only filled for structured
    // bricks
    inline const brickTopologyMap& map() const { return map_; }
};

}

}

#endif
