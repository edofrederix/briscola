#ifndef uniformMesh_H
#define uniformMesh_H

#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

class uniformMesh
:
    public rectilinearMesh
{
    // Uniform cell size
    vector cellSize_;

    // Set uniform metrics
    void setMetrics();

public:

    // Runtime type information
    TypeName("uniform");

    // Construct from IOdictionary
    uniformMesh(const IOdictionary&);

    // Transfer constructor, destroying argument
    uniformMesh(autoPtr<mesh>&);

    // Copy constructor
    uniformMesh(const uniformMesh& msh);

    // Transfer or copy constructor
    uniformMesh(uniformMesh& msh, bool reuse);

    // Clone mesh
    virtual autoPtr<mesh> clone() const
    {
        return autoPtr<mesh>(new uniformMesh(*this));
    }

    // Destructor
    ~uniformMesh();

    // Return the cell size vector
    inline vector cellSize() const
    {
        return cellSize_;
    }

    // Find cells efficiently using algebraic relation
    virtual labelVector findCell(const vector& p, const label l = 0) const;

    //- Disallow default bitwise assignment
    void operator=(const uniformMesh&);
};

}

}

#endif
