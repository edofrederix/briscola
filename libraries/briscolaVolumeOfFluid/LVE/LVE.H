#ifndef LVE_H
#define LVE_H

#include "truncatedHex.H"
#include "truncatedPiped.H"
#include "SortableList.H"
#include "constants.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class vof;

class LVE
{
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
    inline scalar solveBracketedLVE
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

    // Solvers

    inline scalar exactCubicSolver
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

    inline scalar newtonCubicSolver
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

    inline scalar bisectionCubicSolver
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

public:

    // Construct from rectilinear bool
    LVE(const bool rectilinear);

    // Destructor
    ~LVE();

    // Solve the local volume enforcement problem on arbitrary quadrilater
    // hexahedrons
    inline scalar aLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& f
    ) const;

    // CIBRAVE algorithm

    inline scalar cLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& f,
        const scalar& Vtotal
    ) const;

    inline scalar cLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& f
    ) const;

    // Solve the local volume enforcement problem on parallelepiped quadrilater
    // hexahedrons
    inline scalar pLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& f
    ) const;

    // Solve the local volume enforcement problem
    inline scalar operator()
    (
        const scalar f,
        const vertexVector& v,
        const scalar cv,
        const vector& n
    ) const
    {
        const vector m(n/Foam::mag(n));

        return rectilinear_ ? pLVE(v,m,f) : cLVE(v,m,f,cv);
    }
};

}

}

}

#include "LVEI.H"

#endif
