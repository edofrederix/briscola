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
    // Reference to the Vof object
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

    // Compute only one component of a vectorial product
    inline scalar vectorialProductProjection
    (
        const vector& v1,
        const vector& v2,
        const int& i
    ) const
    {
        if (i == 0)
        {
            return v1[1] * v2[2] - v1[2] * v2[1];
        }
        else if (i == 1)
        {
            return v1[2] * v2[0] - v1[0] * v2[2];
        }
        else
        {
            return v1[0] * v2[1] - v1[1] * v2[0];
        }
    }

public:

    // Construct from vof object
    LVE(const vof&);

    // Destructor
    ~LVE();

    // Solve the local volume enforcement problem
    scalar operator()(const labelVector&, const vector& n) const;

    // Solve the local volume enforcement problem
    inline scalar operator()
    (
        const label& i,
        const label& j,
        const label& k,
        const vector& n
    ) const
    {
        return this->operator()(labelVector(i,j,k), n);
    }

    // Solve the local volume enforcement problem on arbitrary quadrilater
    // hexahedrons
    scalar aLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& f
    ) const;

    // CIBRAVE algorithm

    scalar cLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& f,
        const scalar& Vtotal
    ) const;

    scalar cLVE
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

    scalar exactCubicSolver
    (
        const scalar& Vf,
        const scalar& alpha0,
        const scalar& alpha1,
        const scalar& alpha2,
        const scalar& alpha3,
        const scalar& Cmin,
        const scalar& Cmax,
        const scalar& Vkmin,
        const scalar& Vkmax,
        const scalar& Vtotal
    ) const;

    scalar newtonCubicSolver
    (
        const scalar& Vf,
        const scalar& alpha0,
        const scalar& alpha1,
        const scalar& alpha2,
        const scalar& alpha3,
        const scalar& Cmin,
        const scalar& Cmax,
        const scalar& Vkmin,
        const scalar& Vkmax,
        const scalar& Vtotal
    ) const;

    scalar bisectionCubicSolver
    (
        const scalar& Vf,
        const scalar& alpha0,
        const scalar& alpha1,
        const scalar& alpha2,
        const scalar& alpha3,
        const scalar& Cmin,
        const scalar& Cmax,
        const scalar& Vkmin,
        const scalar& Vkmax,
        const scalar& Vtotal
    ) const;

};

}

}

}

#endif
