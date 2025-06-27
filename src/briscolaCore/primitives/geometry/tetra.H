#ifndef tetra_H
#define tetra_H

#include "vector.H"
#include "FixedList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "contiguous.H"
#include "utils.H"

namespace Foam
{

namespace briscola
{

// Forward declarations
class tessellation;

/**

@brief Tetrahedral geometric object

Simplest three-dimensional object with four vertices and four faces each with
three vertices. Some basic functions such as volume computation and truncation
are defined for tets. Using tet tessellations, one can compute such basic
functions also for more complex geometrical objects. The class is called 'tetra'
and not 'tet' to be consistent with the 'hexa' class which cannot be called
'hex', see hexa.H.

*/

class tetra
{
protected:

    // Vertices
    FixedList<vector,4> v_;

    // Edge interpolation function needed to compute truncations
    inline vector interp
    (
        const scalarList& d,
        const label i,
        const label j
    ) const
    {
        scalar t = d[i]/(d[i] - d[j]);
        return v_[i] + t*(v_[j] - v_[i]);
    }

public:

    // Null construct
    inline tetra()
    {}

    // Copy constructor
    inline tetra(const tetra& t)
    :
        v_(t.v_)
    {}

    // Construct from four vertices
    inline tetra
    (
        const vector& lba,
        const vector& lta,
        const vector& rba,
        const vector& lbf
    )
    {
        v_[0] = lba;
        v_[1] = lta;
        v_[2] = rba;
        v_[3] = lbf;

        #ifdef FULLDEBUG
        if (volume() < -1e-12)
            FatalErrorInFunction
                << "Tet has negative volume" << endl << abort(FatalError);
        #endif
    }

    // Construct from a list of four vertices
    inline tetra(const vectorList& v)
    :
        v_(v)
    {
        #ifdef FULLDEBUG

        if (volume() < -1e-12)
            FatalErrorInFunction
                << "Tet has negative volume" << endl << abort(FatalError);

        #endif
    }

    // Construct from a fixed list of four vertices
    inline tetra(const FixedList<vector,4>& v)
    :
        v_(v)
    {
        #ifdef FULLDEBUG

        if (volume() < -1e-12)
            FatalErrorInFunction
                << "Tet has negative volume" << endl << abort(FatalError);

        #endif
    }

    // Access to the vertices
    inline FixedList<vector,4>& v() { return v_; }
    inline const FixedList<vector,4>& v() const { return v_; }

    // Named vertices

    inline vector& lba() { return v_[0]; }
    inline vector& lta() { return v_[1]; }
    inline vector& rba() { return v_[2]; }
    inline vector& lbf() { return v_[3]; }

    inline const vector& lba() const { return v_[0]; }
    inline const vector& lta() const { return v_[1]; }
    inline const vector& rba() const { return v_[2]; }
    inline const vector& lbf() const { return v_[3]; }

    // Compute the center as the average of the vertices
    inline vector center() const
    {
        vector c = v_[0]/4.0;

        for (label i = 1; i < 4; i++)
            c += v_[i]/4.0;

        return c;
    }

    // Compute the face center as the average of the four face vertices
    inline vector faceCenter(const label i) const
    {
        #ifdef FULLDEBUG
        if (i < 0 || i > 3)
            FatalErrorInFunction
                << "Invalid face number" << endl << abort(FatalError);
        #endif

        return
            v_[vertexNumsInTetFace[i][0]]/3.0
          + v_[vertexNumsInTetFace[i][1]]/3.0
          + v_[vertexNumsInTetFace[i][2]]/3.0;
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
            v_[vertexNumsInTetEdge[i][0]]/2.0
          + v_[vertexNumsInTetEdge[i][1]]/2.0;
    }

    // Face area normals, outward pointing
    inline vector faceAreaNormal(const label i) const
    {
        #ifdef FULLDEBUG
        if (i < 0 || i > 3)
            FatalErrorInFunction
                << "Invalid face number" << endl << abort(FatalError);
        #endif

        return triAreaVector
        (
            v_[vertexNumsInTetFace[i][0]],
            v_[vertexNumsInTetFace[i][1]],
            v_[vertexNumsInTetFace[i][2]]
        );
    }

    // Face area
    inline scalar faceArea(const label i) const
    {
        return Foam::mag(faceAreaNormal(i));
    }

    // Compute the volume using the tetVolume function
    inline scalar volume() const
    {
        return tetVolume(v_[0], v_[1], v_[2], v_[3]);
    }

    // Compute truncated tessellation from this tet and a plane defined by
    // dot(x,n) + C = 0. By convention, the normal n points away from the
    // resulting truncated volume.
    tessellation truncate(const vector& n, const scalar C) const;

    // Compute intersection with another object
    template<class Type>
    tessellation intersect(const Type& t) const;

    // Compute intersection between this tet and another tet, by truncating it
    // over all four faces of the given tet
    tessellation intersect(const tetra& t) const;

    // Compute intersection between this tet and a tessellation
    tessellation intersect(const tessellation& T) const;

    // Volume of this tet truncated by a geometrical plane
    scalar truncationVolume(const vector& n, const scalar C) const;

    // Volume of the intersection with another geometrical object
    template<class Type>
    scalar intersectionVolume(const Type& t) const;
};

}

// Tetra class is contiguous
template<>
inline bool contiguous<briscola::tetra>() { return true; }

}

#endif
