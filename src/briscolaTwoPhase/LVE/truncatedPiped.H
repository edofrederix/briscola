#ifndef truncatedPiped_H
#define truncatedPiped_H

#include "fvMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class truncatedPiped
{
    // Vertices
    vector lba_;
    vector rba_;
    vector lta_;
    vector lbf_;

    // Plane normal
    vector n_;

    // Plane constant
    scalar C_;

    // Disallow default bitwise assignment
    void operator=(const truncatedPiped&);

public:

    // Construct from components
    inline truncatedPiped
    (
        const vector& lba,
        const vector& rba,
        const vector& lta,
        const vector& lbf,
        const vector& n,
        const scalar& C
    );

    // Construct from hexahedron
    inline truncatedPiped
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& C
    );

    // Copy constructor
    inline truncatedPiped(const truncatedPiped&);

    // Clone
    inline autoPtr<truncatedPiped> clone() const
    {
        return autoPtr<truncatedPiped>(new truncatedPiped(*this));
    }

    // Destructor
    inline ~truncatedPiped()
    {}

    // Access to the vertices
    inline vector lba() const { return lba_; }
    inline vector rba() const { return rba_; }
    inline vector lta() const { return lta_; }
    inline vector lbf() const { return lbf_; }

    // Return the plane normal
    inline vector n() const { return n_; }

    // Return the plane constant
    inline scalar C() const { return C_; }

    // Compute the volume
    inline scalar volume() const;

};

}

}

}

#include "truncatedPipedI.H"

#endif
