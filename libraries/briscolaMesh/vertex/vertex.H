#ifndef vertex_H
#define vertex_H

#include "meshObject.H"
#include "PtrList.H"
#include "vector.H"

namespace Foam
{

namespace briscola
{

class edge;

class vertex
:
    public meshObject<edge>,
    public vector
{
    // Index of this vertex in the geometry data
    const label j_;

public:

    // Construct from components
    vertex(const edge&, const label, const label);

    // Copy constructor
    vertex(const vertex&);

    // Clone this vertex (needed by PtrList)
    inline autoPtr<vertex> clone() const
    {
        return autoPtr<vertex>(new vertex(*this));
    }

    // Destructor
    virtual ~vertex();

    // Access to the edge
    inline const edge& parentEdge() const
    {
        return parent_;
    }

    // Access to the index of this vertex in the geometry data
    inline label j() const
    {
        return j_;
    }

    // Access to the vertex
    inline vector& v()
    {
        return *this;
    }

    // Access to the vertex
    inline const vector& v() const
    {
        return *this;
    }

    // Assignment operator
    void operator=(const vertex&);

    //- Ostream Operator
    friend Ostream& operator<<(Ostream&, const vertex&);
};

}

}

#endif
