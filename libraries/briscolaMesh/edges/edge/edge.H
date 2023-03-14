#ifndef edge_H
#define edge_H

#include "meshObject.H"
#include "PtrList.H"
#include "block.H"
#include "vertex.H"
#include "scalarField.H"
#include "vectorField.H"

namespace Foam
{

namespace briscola
{

class edge
:
    public meshObject<face>
{
    // Block of vertices
    labelBlock v_;

    // Number of cells on this edge
    labelVector N_;

    // List of vertices
    PtrList<vertex> vertices_;

    // Disallow default bitwise assignment
    void operator=(const edge&);

public:

    // Runtime type information
    TypeName("edge");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        edge,
        dictionary,
        (
            const face& f,
            const label num,
            const labelBlock& vertexNums,
            const labelVector& N
        ),
        (f, num, vertexNums, N)
    );

    // Construct from components
    edge(const face&, const label, const labelBlock&, const labelVector&);

    // Copy constructor
    edge(const edge&);

    // Copy constructor with reference to new face
    edge(const edge&, const face& f);

    // Clone edge with reference to new face
    virtual autoPtr<edge> clone(const face&) const = 0;

    // Selector
    static autoPtr<edge> New
    (
        const face&,
        const label,
        const labelBlock&,
        const labelVector&
    );

    // Destructor
    virtual ~edge();

    // Access to the face
    inline const face& parentFace() const
    {
        return parent_;
    }

    // Access to the list of vertices
    inline const PtrList<vertex>& vertices() const
    {
        return vertices_;
    }

    // Number of cells on this edge
    inline labelVector N() const
    {
        return N_;
    }

    // Zeroth vertex
    inline const vertex& v0() const
    {
        return vertices_[0];
    }

    // First vertex
    inline const vertex& v1() const
    {
        return vertices_[1];
    }

    // Block of vertices
    inline const labelBlock& v() const
    {
        return v_;
    }

    // Distance vector of the edge
    inline vector dist() const
    {
        return v1()-v0();
    }

    // Length of the edge
    virtual scalar length() const = 0;

    // Parametrized coordinates
    virtual tmp<vectorField> operator()(const scalarField&) const = 0;

    // Dictionary of this edge
    dictionary dict() const;

    // Equality operator
    bool operator==(const edge&) const;

    // Unequality operator
    bool operator!=(const edge&) const;

    //- Ostream Operator
    friend Ostream& operator<<(Ostream&, const edge&);
};

}

}

#endif
