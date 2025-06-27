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

/**

@brief Local Volume Enforcement (LVE) solver

This class offers tools to solve the LVE problem. The LVE problem is defined by
a cell, a plane normal and a volume fraction which is enclosed by the plane and
the faces of the cell. Given the plane normal, the LVE's task is to compute the
scalar parameter in the equations dot(x,n) + C = 0 that defines the plane, with
normal n. For parallelepiped quadrilateral hexahedrons an analytical approach
can be used. For general hexahedrons the CIBRAVE algorithm of López et al.
(2016) is used.

*/

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

    // Solvers

    inline scalar exactCubicSolver
    (
        const scalar Vf,
        const scalarList& alpha,
        const scalar Cmin,
        const scalar Cmax,
        const scalar Wmin,
        const scalar Wmax,
        const scalar V
    ) const;

    inline scalar NewtonCubicSolver
    (
        const scalar Vf,
        const scalarList& alpha,
        const scalar Cmin,
        const scalar Cmax,
        const scalar Wmin,
        const scalar Wmax,
        const scalar V
    ) const;

    inline scalar bisectionCubicSolver
    (
        const scalar Vf,
        const scalarList& alpha,
        const scalar Cmin,
        const scalar Cmax,
        const scalar Wmin,
        const scalar Wmax,
        const scalar V
    ) const;

    // CIBRAVE algorithm

    inline scalar cLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar f,
        const scalar V,
        const vectorList& triAreas,
        const vectorList& triNormals,
        const scalarList& Cj,
        const labelList& m
    ) const;

public:

    // Construct from rectilinear bool
    inline LVE(const bool rectilinear)
    :
        rectilinear_(rectilinear)
    {}

    // Destructor
    inline ~LVE()
    {}

    // CIBRAVE algorithm

    inline scalar cLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar f
    ) const;

    // Solve the local volume enforcement problem on parallelepiped
    // quadrilateral hexahedrons
    inline scalar pLVE
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar f
    ) const;

    // Solve the local volume enforcement problem
    inline scalar operator()
    (
        const scalar f,
        const vertexVector& v,
        const vector& n
    ) const
    {
        const vector m(n/Foam::mag(n));

        return rectilinear_ ? pLVE(v,m,f) : cLVE(v,m,f);
    }

    scalar fluxVolumeLVE
    (
        const vertexVector& v,
        const scalar V,
        const scalar Vf,
        const label u,
        const label l
    ) const;
};

}

}

}

#include "LVEI.H"

#endif
