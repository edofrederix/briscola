#ifndef immersedBoundary_H
#define immersedBoundary_H

#include "immersedBoundaryFwd.H"

#include "List.H"
#include "immersedBoundaryShape.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "fvMeshMetrics.H"
#include "meshField.H"
#include "linearSystem.H"
#include "rectilinearMesh.H"
#include "utils.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class immersedBoundary
{
protected:

    // Reference to mesh
    const fvMesh& fvMsh_;

    // Reference to fvMeshMetrics
    const fvMeshMetrics<MeshType>& fvMshMetrics_;

    // Name of the immersed boundary
    const word name_;

    // Switch for whether there are overlapping shapes within the mesh
    Switch shapeOverlap_;

    // The IB mask is equal to 0 in cell centers outside the IB, and equal to 1
    // in cell centers inside the IB.
    meshField<label,MeshType> mask_;

    // The ghost cell mask is equal to 1 in IB cells with at least one
    // neighboring fluid cell and 0 in other cells
    meshField<label,MeshType> ghostMask_;

    // The wall adjacent mask is equal to 1 in fluid cells with at least one
    // neighboring cell in the IB, and 0 in other cells
    meshField<label,MeshType> wallAdjMask_;

    // This wall distance field stores the wall distances from IB-adjacent cells
    // to the IB as faceScalars
    faceField<scalar,MeshType> wallDistAdj_;

    // This wall distance field stores the wall distances from ghost cells to
    // the IB as faceScalars
    faceField<scalar,MeshType> wallDistGhost_;

    // The neighbor distance field stores the distances to the second neighbor
    // fluid cell for the quadratic fit on non-uniform grids
    faceField<scalar,MeshType> neighborDist_;

    // Field containing the coordinates of the ghost mirror points stored in the
    // ghost cells
    meshField<vector,MeshType> mirrorPoints_;

    // List of shapes making up the immersed boundary
    FastPtrList<immersedBoundaryShape> shapes_;

    // Set mask fields
    void setMasks();

    // Calculate wall distance fields
    void calculateWallDistances();

    // Set mirror points
    void setMirrorPoints();

public:

    // Declare name of the class and its debug switch
    TypeName("immersedBoundary");

    // Construct from fvMeshMetrics and solverDict
    immersedBoundary(const fvMeshMetrics<MeshType>&, const dictionary&);

    // Destructor
    virtual ~immersedBoundary();

    // Returns true if the position vector (x,y,z) is inside the immersed
    // boundary
    bool isInside(vector) const;

    // Returns dimensionless wall distance given two points The first point is a
    // fluid cell and the second point is an IB cell
    scalar wallDistance(vector, vector) const;

    // Returns the distance from a ghost cell to the immersed boundary in the
    // wall-normal direction
    scalar wallNormalDistance(vector) const;

    // Returns the coordinates of the mirror point in the wall-normal direction
    // of a point located inside of the immersed boundary
    vector mirrorPoint(vector) const;

    // Check if any of the IB shapes overlap within the mesh
    void checkOverlap();

    // Inline access functions
    inline const word name() const { return name_; }
    inline Switch shapeOverlap() const { return shapeOverlap_; }

    // Access to the immersed boundary fields

    inline const meshField<label,MeshType>& mask() const
    {
        return mask_;
    }

    inline const meshField<label,MeshType>& ghostMask() const
    {
        return ghostMask_;
    }

    inline const meshField<label,MeshType>& wallAdjMask() const
    {
        return wallAdjMask_;
    }

    inline const faceField<scalar,MeshType>& wallDistAdj() const
    {
        return wallDistAdj_;
    }

    inline const faceField<scalar,MeshType>& wallDistGhost() const
    {
        return wallDistGhost_;
    }

    inline const faceField<scalar,MeshType>& neighborDist() const
    {
        return neighborDist_;
    }

    inline const meshField<vector,MeshType>& mirrorPoints() const
    {
        return mirrorPoints_;
    }
};

}

}

}

#endif
