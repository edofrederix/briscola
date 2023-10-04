#ifndef truncatedTet_H
#define truncatedTet_H

#include "fvMesh.H"
#include "vectorList.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class truncatedTet
{
    // Vertices
    vectorList v_;

    // Plane normal
    vector n_;

    // Plane constant
    scalar C_;

    // List of vertex constants
    scalarList phi_;

    // Set the vertex constants
    inline void setVertexConstants();

    // Function to compute the intersection of the plane and the line segment
    // between two vertices
    inline vector intersectionPoint(const label i, const label j) const
    {
        const vector d = v_[j] - v_[i];
        scalar factor = - (C_ + (n_ & v_[i]))/stabilise(n_ & d, 1e-16);

        if (factor >= 0 && factor <= 1)
            return v_[i] + factor * d;
        else if (factor < 0)
            factor = 1e-12;
        else if (factor > 1)
            factor = 1 - 1e-12;
        else
            factor = 0.5;

        return v_[i] + factor * d;
    }

    // Disallow default bitwise assignment
    void operator=(const truncatedTet&);

public:

    // Construct from components
    inline truncatedTet
    (
        const vectorList& v,
        const vector& n,
        const scalar& C
    );

    // Copy constructor
    inline truncatedTet(const truncatedTet&);

    // Clone
    inline autoPtr<truncatedTet> clone() const
    {
        return autoPtr<truncatedTet>(new truncatedTet(*this));
    }

    // Destructor
    virtual ~truncatedTet();

    // Access to the vertices
    inline const vectorList& v() const { return v_; }

    // Return the plane normal
    inline vector n() const { return n_; }

    // Return the plane constant
    inline scalar C() const { return C_; }

    // Access to the vertex constants
    inline const scalarList& phi() const { return phi_; }

    // Compute the volume
    inline scalar volume() const;
};

}

}

}

#include "truncatedTetI.H"

#endif
