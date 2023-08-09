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
    vector v_lba_;
    vector v_rba_;
    vector v_lta_;
    vector v_lbf_;

    // Plane normal
    vector n_;

    // Plane constant
    scalar C_;

    // Disallow default bitwise assignment
    void operator=(const truncatedPiped&);

public:

    // Construct from components
    truncatedPiped
    (
        const vector& v_lba,
        const vector& v_rba,
        const vector& v_lta,
        const vector& v_lbf,
        const vector& n,
        const scalar& C
    );

    // Construct from hexahedron
    truncatedPiped
    (
        const vertexVector& vertices,
        const vector& n,
        const scalar& C
    );

    // Copy constructor
    truncatedPiped(const truncatedPiped&);

    // Clone
    autoPtr<truncatedPiped> clone() const
    {
        return autoPtr<truncatedPiped>(new truncatedPiped(*this));
    }

    // Destructor
    virtual ~truncatedPiped();

    // Access to the vertices
    inline vector v_lba() const { return v_lba_; }
    inline vector v_rba() const { return v_rba_; }
    inline vector v_lta() const { return v_lta_; }
    inline vector v_lbf() const { return v_lbf_; }

    // Return the plane normal
    inline vector n() const { return n_; }

    // Return the plane constant
    inline scalar C() const { return C_; }

    // Compute the volume
    virtual scalar volume() const;

};

}

}

}

#endif
