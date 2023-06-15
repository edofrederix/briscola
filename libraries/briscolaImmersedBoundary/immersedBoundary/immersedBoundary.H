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

class immersedBoundary
{
    // Mesh
    const fvMesh& fvMsh_;

    // The colocated IB mask is equal to 0 in colocated cell
    // centers outside the IB, and equal to 1 in colocated cell
    // centers inside the IB.
    colocatedScalarField colMask_;

    // The staggered IB mask is equal to 0 in staggered cell
    // centers outside the IB, and equal to 1 in staggered cell
    // centers inside the IB.
    staggeredScalarField stagMask_;

    // List of shapes making up the immersed boundary
    PtrList<shape> shapes_;

public:

    // Construct from solverDict and fvMesh
    immersedBoundary
    (
        IOdictionary,
        const fvMesh&
    );

    // Destructor
    virtual ~immersedBoundary();

    // Returns true if the position vector (x,y,z)
    // is inside the immersed boundary
    bool isInside(vector);

    // Return functions for mask fields
    inline colocatedScalarField& colMask() { return colMask_; }
    inline staggeredScalarField& stagMask() { return stagMask_; }

};

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam

#endif
