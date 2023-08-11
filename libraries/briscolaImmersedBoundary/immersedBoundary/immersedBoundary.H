#ifndef immersedBoundary_H
#define immersedBoundary_H

#include "List.H"
#include "cylinder.H"
#include "sphere.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "fvMeshMetrics.H"
#include "meshField.H"

namespace Foam
{

namespace briscola
{

namespace ibm
{

using namespace fv;

template<class MeshType>
class immersedBoundary
{
    // Mesh
    const fvMesh& fvMsh_;

    // Stability factor for xi. If xi > stability factor,
    // the IBM will revert to a more stable first order fit.
    scalar xiStabilityFactor_;

    // The IB mask is equal to 0 in cell centers outside the IB,
    // and equal to 1 in cell centers inside the IB.
    meshField<scalar,MeshType> mask_;

    // The wall adjacent mask is equal to 1 in fluid cells with at
    // least one neighboring cell in the IB, and 0 in other cells
    meshField<scalar,MeshType> wallAdjMask_;

    // The wall distance field stores the wall distances to the
    // IB in the components of the stencils
    meshField<stencil,MeshType> wallDist_;

    // List of shapes making up the immersed boundary
    PtrList<shape> shapes_;

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

    // Return functions for mask fields
    inline meshField<scalar,MeshType>& mask() { return mask_; }
    inline meshField<scalar,MeshType>& wallAdjMask() { return wallAdjMask_; }
    inline meshField<stencil,MeshType>& wallDist() { return wallDist_; }

    // Return stability factor
    inline scalar xiStabilityFactor() { return xiStabilityFactor_; }

};

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam

#endif
