#ifndef immersedBoundary_H
#define immersedBoundary_H

#include "immersedBoundaryFwd.H"

#include "List.H"
#include "shape.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "fvMeshMetrics.H"
#include "meshField.H"
#include "linearSystem.H"
#include "rectilinearMesh.H"
#include "blockUtils.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class immersedBoundary
{
    // Reference to mesh
    const fvMesh& fvMsh_;

    // IBM type
    word type_;

    // Mass source active switch
    bool massSourceActive_;

    // The IB mask is equal to 0 in cell centers outside the IB,
    // and equal to 1 in cell centers inside the IB.
    meshField<label,MeshType> mask_;

    // The wall adjacent mask is equal to 1 in fluid cells with at
    // least one neighboring cell in the IB, and 0 in other cells
    meshField<label,MeshType> wallAdjMask_;

    // The ghost cell mask is equal to 1 in IB cells with at least
    // one neighboring fluid cell and 0 in other cells
    meshField<label,MeshType> ghostMask_;

    // The wall distance field stores the wall distances to the
    // IB in the components of the stencils
    meshField<faceScalar,MeshType> wallDist_;

    // The neighbor distance field stores the distances to neighboring
    // fluid cell for the quadratic fit on non-uniform grids
    meshField<faceScalar,MeshType> neighborDist_;

    // List of shapes making up the immersed boundary
    PtrList<shape> shapes_;

    // List of cell coordinates to be sent to/received from
    // processors.
    // sendRecvList[pr][ps] contains the list of cell indices
    // that processor pr needs to receive from processor ps
    // List<List<List<List<List<labelVector>>>>> sendRecvList_;

public:

    // Declare name of the class and its debug switch
    ClassName("immersedBoundary");

    // Construct from solverDict and fvMesh
    immersedBoundary
    (
        IOdictionary,
        const fvMesh&
    );

    // Destructor
    ~immersedBoundary();

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

    // Apply penalization method to a linear system
    void penalization(linearSystem<stencil,Type,MeshType>&);

    // Apply Fadlun immersed boundary method for Dirichlet BC to a linear system
    void FadlunIBM(linearSystem<stencil,Type,MeshType>&);

    // Return colocated volume source term for a staggered velocity field.
    // Only defined for staggeredScalarImmersedBoundary
    tmp<meshField<scalar,colocated>> IBMSource
    (
        const meshField<scalar,staggered>&
    );

    // Inline access functions
    inline word& type() { return type_; }
    inline bool massSourceActive() {return massSourceActive_; }
    inline meshField<label,MeshType>& mask() { return mask_; }
    inline meshField<label,MeshType>& wallAdjMask() { return wallAdjMask_; }
    inline meshField<label,MeshType>& ghostMask() { return ghostMask_; }
    inline meshField<faceScalar,MeshType>& wallDist() { return wallDist_; }
    inline meshField<faceScalar,MeshType>& neighborDist() { return neighborDist_; }

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
