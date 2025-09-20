#ifndef face_H
#define face_H

#include "meshObject.H"
#include "PtrList.H"
#include "edge.H"
#include "block.H"

namespace Foam
{

namespace briscola
{

class brick;

class face
:
    public meshObject<brick>
{
    // Block of vertices
    labelBlock v_;

    // Number of cells on this face
    labelVector N_;

    // List of edges
    PtrList<edge> edges_;

    // Create edges
    void createEdges();

public:

    // Construct from components
    face(const brick&, const label, const labelBlock&, const labelVector&);

    // Copy constructor
    face(const face&);

    // Copy constructor with reference to new brick
    face(const face&, const brick&);

    // Clone face with reference to new brick
    inline autoPtr<face> clone(const brick& b) const
    {
        return autoPtr<face>(new face(*this, b));
    }

    // Destructor
    virtual ~face();

    // Access to the brick
    inline const brick& parentBrick() const
    {
        return parent_;
    }

    // Access to the list of edges
    inline const PtrList<edge>& edges() const
    {
        return edges_;
    }

    // Number of cells on this face
    inline labelVector N() const
    {
        return N_;
    }

    // Zeroth edge
    inline const edge& e0() const
    {
        return edges_[0];
    }

    // First edge
    inline const edge& e1() const
    {
        return edges_[1];
    }

    // Second edge
    inline const edge& e2() const
    {
        return edges_[2];
    }

    // Third edge
    inline const edge& e3() const
    {
        return edges_[3];
    }

    // Zeroth vertex
    inline const vertex& v0() const
    {
        return e0().v0();
    }

    // First vertex
    inline const vertex& v1() const
    {
        return e0().v1();
    }

    // Second vertex
    inline const vertex& v2() const
    {
        return e1().v0();
    }

    // Third vertex
    inline const vertex& v3() const
    {
        return e1().v1();
    }

    // Block of vertices
    inline const labelBlock& v() const
    {
        return v_;
    }

    // Area of the face (using area of triangles)
    inline scalar area() const
    {
        return
            triArea(v0(), v1(), v2())
          + triArea(v3(), v2(), v1());
    }

    // Return list of vertices
    inline labelList vertexList() const
    {
        labelList vertices;
        forAllBlockLinear(v_, i)
            vertices.append(v_(i));

        return vertices;
    }

    // Equality operator
    bool operator==(const face&) const;

    // Unequality operator
    bool operator!=(const face&) const;

    //- Ostream Operator
    friend Ostream& operator<<(Ostream&, const face&);
};

}

}

#endif
