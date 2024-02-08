#ifndef immersedBoundary_H
#define immersedBoundary_H

#include "immersedBoundaryFwd.H"

#include "List.H"
#include "shape.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "fvMeshMetrics.H"
#include "meshField.H"
#include "linearSystem.H"
#include "rectilinearMesh.H"
#include "briscolaUtils.H"

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

    // The IB mask is equal to 0 in cell centers outside the IB,
    // and equal to 1 in cell centers inside the IB.
    meshField<label,MeshType> mask_;

    // The ghost cell mask is equal to 1 in IB cells with at least
    // one neighboring fluid cell and 0 in other cells
    meshField<label,MeshType> ghostMask_;

    // The wall adjacent mask is equal to 1 in fluid cells with at
    // least one neighboring cell in the IB, and 0 in other cells
    meshField<label,MeshType> wallAdjMask_;

    // List of shapes making up the immersed boundary
    PtrList<shape> shapes_;

public:

    // Declare name of the class and its debug switch
    TypeName("immersedBoundary");

    // Construct from solverDict and fvMesh
    immersedBoundary(const fvMesh&);

    // Destructor
    virtual ~immersedBoundary();

    // Returns true if the position vector (x,y,z)
    // is inside the immersed boundary
    bool isInside(vector);

    // Returns dimensionless wall distance given two points
    // The first point is a fluid cell and the second point is an IB cell
    scalar wallDistance(vector, vector);

    // Returns the distance from a ghost cell to the immersed boundary in
    // the wall-normal direction
    scalar wallNormalDistance(vector);

    // Returns the coordinates of the mirror point in the wall-normal direction
    // of a point located inside of the immersed boundary
    vector mirrorPoint(vector);

    // Inline access functions
    inline const meshField<label,MeshType>& mask() const { return mask_; }
    inline const meshField<label,MeshType>& ghostMask() const { return ghostMask_; }
    inline const meshField<label,MeshType>& wallAdjMask() const { return wallAdjMask_; }

};

}

}

}

#endif
