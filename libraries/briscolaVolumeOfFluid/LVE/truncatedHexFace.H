#ifndef truncatedHexFace_H
#define truncatedHexFace_H

#include "vertexVector.H"
#include "vectorList.H"


namespace Foam
{

namespace briscola
{

namespace fv
{

class truncatedHexFace
{
protected:

    // List of vertices
    vectorList v_;

public:

    // Construct from components
    truncatedHexFace(const vectorList&);

    // Null constructor
    truncatedHexFace()
    {}

    // Copy constructor
    truncatedHexFace(const truncatedHexFace&);

    // Assignment operator
    void operator=(const truncatedHexFace&);

    // Clone
    autoPtr<truncatedHexFace> clone() const
    {
        return autoPtr<truncatedHexFace>(new truncatedHexFace(*this));
    }

    // Destructor
    virtual ~truncatedHexFace();

    // Compute vector area as an average of all the triangular decompositions of
    // the plane
    virtual vector vectorArea() const;

    // Compute outward normal
    virtual vector normal() const;

    // Compute the area
    virtual scalar area() const;

    // Return the vertices
    inline const vectorList& v() const
    {
        return v_;
    }
};

}

}

}

#endif
