#ifndef LVE_H
#define LVE_H

#include "fvMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class vof;

class LVE
{
    // Reference to the VoF object
    const vof& vf_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Rectilinear bool (stored here for efficiency)
    bool rectilinear_;

    // Tolerance
    static const scalar tol_;

    // Max number of iterations
    static const label maxIter_;

    // Disallow default bitwise assignment
    void operator=(const LVE&);

    // Disallow copy constructor
    LVE(const LVE&);

    // Functions to solve the local volume enforcement problem given a solution
    // bracket
    scalar solveBracketedLVE
    (
        const vertexVector& v,
        const scalar& V,
        const vector& n,
        const scalar& f,
        scalar CMin,
        scalar CMax,
        scalar fMin,
        scalar fMax
    ) const;

public:

    // Construct from vof object
    LVE(const vof&);

    // Destructor
    ~LVE();

    // Solve the local volume enforcement problem
    scalar operator()(const labelVector&, const vector& n) const;

    // Solve the local volume enforcement problem on arbitrary quadrilater
    // hexahedrons
    scalar aLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& f
    ) const;

    // Solve the local volume enforcement problem on parallelepiped quadrilater
    // hexahedrons
    scalar pLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& f
    ) const;
};

}

}

}

#endif
