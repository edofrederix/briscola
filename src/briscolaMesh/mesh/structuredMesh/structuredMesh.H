#ifndef structuredMesh_H
#define structuredMesh_H

#include "unstructuredMesh.H"

namespace Foam
{

namespace briscola
{

class structuredMesh
:
    public unstructuredMesh
{
    // Global mesh size
    labelVector N_;

    // Set global mesh size
    void setGlobalMeshSize();

public:

    // Runtime type information
    TypeName("structured");

    // Construct from IOdictionary
    structuredMesh(const IOdictionary&);

    // Transfer constructor, destroying argument
    structuredMesh(autoPtr<mesh>&);

    // Copy constructor
    structuredMesh(const structuredMesh& msh);

    // Clone mesh
    virtual autoPtr<mesh> clone() const
    {
        return autoPtr<mesh>(new structuredMesh(*this));
    }

    // Destructor
    ~structuredMesh();

    // Global mesh size
    inline labelVector N() const
    {
        return N_;
    }

    //- Disallow default bitwise assignment
    void operator=(const structuredMesh&);
};

}

}

#endif
