#ifndef hexa_H
#define hexa_H

#include "vector.H"
#include "vectorList.H"
#include "vertexVector.H"
#include "contiguous.H"

namespace Foam
{

namespace briscola
{

// Forward declarations
class tessellation;
class piped;

/**

@brief Hexahedral geometric object

Very simple hexahedral object that allows to compute geometric properties such
as its center, face centers and volume, and that allows for truncation with
other geometrical objects such as a plane or a tetrahedral. The reason that this
class is called 'hexa' and not 'hex' is because OpenFOAM has already defined the
hex class.

*/

class hexa
{
protected:

    // Vertices
    vertexVector v_;

public:

    // Null construct
    inline hexa()
    {}

    // Copy constructor
    inline hexa(const hexa& h)
    :
        v_(h.v_)
    {}

    // Construct from eight vertices
    inline hexa
    (
        const vector& lba,
        const vector& rba,
        const vector& lta,
        const vector& rta,
        const vector& lbf,
        const vector& rbf,
        const vector& ltf,
        const vector& rtf
    )
    {
        v_.lba() = lba;
        v_.rba() = rba;
        v_.lta() = lta;
        v_.rta() = rta;
        v_.lbf() = lbf;
        v_.rbf() = rbf;
        v_.ltf() = ltf;
        v_.rtf() = rtf;
    }

    // Construct from vertex vector
    inline hexa(const vertexVector& v)
    :
        v_(v)
    {}

    // Construct from piped
    hexa(const piped& p);

    // Construct from a list of eight vertices
    inline hexa(const vectorList& v)
    {
        v_.lba() = v[0];
        v_.rba() = v[1];
        v_.lta() = v[2];
        v_.rta() = v[3];
        v_.lbf() = v[4];
        v_.rbf() = v[5];
        v_.ltf() = v[6];
        v_.rtf() = v[7];
    }

    // Access to the vertices
    inline vertexVector& v() { return v_; }
    inline const vertexVector& v() const { return v_; }

    // Named vertices

    inline vector& lba() { return v_[0]; }
    inline vector& rba() { return v_[1]; }
    inline vector& lta() { return v_[2]; }
    inline vector& rta() { return v_[3]; }
    inline vector& lbf() { return v_[4]; }
    inline vector& rbf() { return v_[5]; }
    inline vector& ltf() { return v_[6]; }
    inline vector& rtf() { return v_[7]; }

    inline const vector& lba() const { return v_[0]; }
    inline const vector& rba() const { return v_[1]; }
    inline const vector& lta() const { return v_[2]; }
    inline const vector& rta() const { return v_[3]; }
    inline const vector& lbf() const { return v_[4]; }
    inline const vector& rbf() const { return v_[5]; }
    inline const vector& ltf() const { return v_[6]; }
    inline const vector& rtf() const { return v_[7]; }

    // Compute the center as the average of the vertices
    inline vector center() const
    {
        vector c = v_[0]/8.0;

        for (label i = 1; i < 8; i++)
            c += v_[i]/8.0;

        return c;
    }

    // Compute the face center as the average of the four face vertices
    inline vector faceCenter(const label i) const
    {
        #ifdef FULLDEBUG
        if (i < 0 || i > 5)
            FatalErrorInFunction
                << "Invalid face number" << endl << abort(FatalError);
        #endif

        return
            v_[vertexNumsInFaceCC[i][0]]/4.0
          + v_[vertexNumsInFaceCC[i][1]]/4.0
          + v_[vertexNumsInFaceCC[i][2]]/4.0
          + v_[vertexNumsInFaceCC[i][3]]/4.0;
    }

    // Compute the edge center as the average of the two edge vertices
    inline vector edgeCenter(const label i) const
    {
        #ifdef FULLDEBUG
        if (i < 0 || i > 11)
            FatalErrorInFunction
                << "Invalid edge number" << endl << abort(FatalError);
        #endif

        return
            v_[vertexNumsInEdge[i][0]]/2.0
          + v_[vertexNumsInEdge[i][1]]/2.0;
    }

    // Face area normals, outward pointing
    inline vector faceAreaNormal(const label i) const
    {
        #ifdef FULLDEBUG
        if (i < 0 || i > 5)
            FatalErrorInFunction
                << "Invalid face number" << endl << abort(FatalError);
        #endif

        return
            0.5
          * (
                (
                    v_[vertexNumsInFaceCC[i][2]]
                  - v_[vertexNumsInFaceCC[i][0]]
                )
              ^ (
                    v_[vertexNumsInFaceCC[i][3]]
                  - v_[vertexNumsInFaceCC[i][1]]
                )
            );
    }

    // Face area
    inline scalar faceArea(const label i) const
    {
        return Foam::mag(faceAreaNormal(i));
    }

    // Compute the volume from the tessellation
    scalar volume() const;

    // Check if hex is valid
    bool valid() const;

    // Compute truncated tessellation from this hex and a plane defined by
    // dot(x,n) + C = 0
    tessellation truncate(const vector& n, const scalar C) const;

    // Compute intersection with another object
    template<class Type>
    tessellation intersect(const Type& t) const;

    // Volume of this hex truncated by a geometrical plane
    scalar truncationVolume(const vector& n, const scalar C) const;

    // Volume of the intersection with another geometrical object
    template<class Type>
    scalar intersectionVolume(const Type& t) const;
};

}

// Hexa class is contiguous
template<>
inline bool contiguous<briscola::hexa>() { return true; }

}

#endif
