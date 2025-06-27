#ifndef piped_H
#define piped_H

#include "FixedList.H"
#include "vector.H"
#include "tensor.H"
#include "vertexVector.H"
#include "contiguous.H"

namespace Foam
{

namespace briscola
{

// Forward declarations
class tessellation;

/**

@brief Parallelepiped geometric object

Very simple parallelepiped object that allows to compute geometric properties
such as its center, face centers and volume, and that allows for truncation with
other geometrical objects such as a plane or a tetrahedral. It is derived from
the hexahedral class, and overwrites a few member functions that can be done
faster without tessellation.

*/

class piped
:
    public hexa
{
protected:

    // Function to check if hex is piped
    inline void checkPiped() const
    {
        if (Foam::mag(rta() - v0() - v1() - p()) > 1e-12)
            FatalErrorInFunction
                << "Invalid piped" << endl << abort(FatalError);

        if (Foam::mag(rbf() - v0() - v2() - p()) > 1e-12)
            FatalErrorInFunction
                << "Invalid piped" << endl << abort(FatalError);

        if (Foam::mag(ltf() - v1() - v2() - p()) > 1e-12)
            FatalErrorInFunction
                << "Invalid piped" << endl << abort(FatalError);

        if (Foam::mag(rtf() - v0() - v1() - v2() - p()) > 1e-12)
            FatalErrorInFunction
                << "Invalid piped" << endl << abort(FatalError);
    }

public:

    // Null construct
    inline piped()
    :
        hexa()
    {}

    // Copy constructor
    inline piped(const piped& p)
    :
        hexa(p)
    {}

    // Construct from base point and orientation tensor
    inline piped(const vector& p, const tensor& T)
    :
        hexa
        (
            p,
            p + T.x(),
            p + T.y(),
            p + T.x() + T.y(),
            p + T.z(),
            p + T.x() + T.z(),
            p + T.y() + T.z(),
            p + T.x() + T.y() + T.z()
        )
    {}

    // Construct from base point and three orientation vectors
    inline piped
    (
        const vector& p,
        const vector& v0,
        const vector& v1,
        const vector& v2
    )
    :
        hexa
        (
            p,
            p + v0,
            p + v1,
            p + v0 + v1,
            p + v2,
            p + v0 + v2,
            p + v1 + v2,
            p + v0 + v1 + v2
        )
    {}

    // Construct from eight vertices
    inline piped
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
    :
        hexa(lba, rba, lta, rta, lbf, rbf, ltf, rtf)
    {
        #ifdef FULLDEBUG
        checkPiped();
        #endif
    }

    // Construct from vertex vector
    inline piped(const vertexVector& v)
    :
        hexa(v)
    {
        #ifdef FULLDEBUG
        checkPiped();
        #endif
    }

    // Construct from a list of eight vertices
    inline piped(const vectorList& v)
    :
        hexa(v)
    {
        #ifdef FULLDEBUG
        checkPiped();
        #endif
    }

    // Base point
    inline vector p() const { return v_[0]; }

    // Orientation vectors
    inline vector v0() const { return v_[1] - p(); }
    inline vector v1() const { return v_[2] - p(); }
    inline vector v2() const { return v_[4] - p(); }

    // Orientation tensor
    inline tensor T() const { return tensor(v0(), v1(), v2()); }

    // Compute the volume without the need for tessellation
    inline scalar volume() const
    {
        return det(T());
    }

    // Volume of this piped truncated by a plane. Analytical function for a
    // piped, without the need for tessellation.
    inline scalar truncationVolume(const vector& m, const scalar C) const;
};

}

// Piped class is contiguous
template<>
inline bool contiguous<briscola::piped>() { return true; }

}

#include "pipedI.H"

#endif
