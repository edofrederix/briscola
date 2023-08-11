#ifndef truncatedHex_H
#define truncatedHex_H

#include "fvMesh.H"
#include "truncatedTet.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class truncatedHex
{
    // Vertices
    vertexVector v_;

    // Plane normal
    vector n_;

    // Plane constant
    scalar C_;

    // List of truncated tets
    PtrList<truncatedTet> tets_;

    // Set the list of truncated tets
    void setTets();

    // Disallow default bitwise assignment
    void operator=(const truncatedHex&);

public:

    // Construct from components
    truncatedHex
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& C
    );

    // Copy constructor
    truncatedHex(const truncatedHex&);

    // Clone
    autoPtr<truncatedHex> clone() const
    {
        return autoPtr<truncatedHex>(new truncatedHex(*this));
    }

    // Destructor
    virtual ~truncatedHex();

    // Access to the vertices
    inline const vertexVector& v() const { return v_; }

    // Return the plane normal
    inline vector n() const { return n_; }

    // Return the plane constant
    inline scalar C() const { return C_; }

    // Access to the faces
    inline const PtrList<truncatedTet>& tets() const { return tets_; }

    // Compute the volume
    virtual scalar volume() const;
};

}

}

}

#endif
