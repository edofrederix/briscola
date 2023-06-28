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

    // The colocated and staggered IB masks are equal to 0 in
    // colocated/staggered cell centers outside the IB, and
    // equal to 1 in cell centers inside the IB.
    colocatedScalarField colMask_;
    staggeredScalarField stagMask_;

    // The colocated and staggered wall adjacent masks are equal
    // to 1 in fluid cells with at least one neighboring cell in
    // the IB, and 0 in other cells
    colocatedScalarField colWallAdjMask_;
    staggeredScalarField stagWallAdjMask_;

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
    inline colocatedScalarField& colWallAdjMask() { return colWallAdjMask_; }
    inline staggeredScalarField& stagWallAdjMask() { return stagWallAdjMask_; }

};

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam

#endif
