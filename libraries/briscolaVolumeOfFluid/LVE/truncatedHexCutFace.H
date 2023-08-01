#ifndef truncatedHexCutFace_H
#define truncatedHexCutFace_H

#include "truncatedHexFace.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class truncatedHexCutFace
:
    public truncatedHexFace
{
public:

    // Construct from components
    truncatedHexCutFace(const vectorList&);

    // Null constructor
    truncatedHexCutFace()
    {}

    // Copy constructor
    truncatedHexCutFace(const truncatedHexCutFace&);

    // Copy constructor
    truncatedHexCutFace(const truncatedHexFace&);

    // Clone
    autoPtr<truncatedHexCutFace> clone() const
    {
        return autoPtr<truncatedHexCutFace>(new truncatedHexCutFace(*this));
    }

    // Destructor
    virtual ~truncatedHexCutFace();

    // Cut face is planar which simplifies the area vector computation
    virtual vector vectorArea() const;
};

}

}

}

#endif
