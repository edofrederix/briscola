#ifndef unstructuredMesh_H
#define unstructuredMesh_H

#include "mesh.H"

namespace Foam
{

namespace briscola
{

class unstructuredMesh
:
    public mesh
{
public:

    // Runtime type information
    TypeName("unstructured");

    // Construct from IOdictionary
    unstructuredMesh(const IOdictionary&);

    // Transfer constructor, destroying argument
    unstructuredMesh(autoPtr<mesh>&);

    // Copy constructor
    unstructuredMesh(const unstructuredMesh& msh);

    // Transfer or copy constructor
    unstructuredMesh(unstructuredMesh& msh, bool reuse);

    // Clone mesh
    virtual autoPtr<mesh> clone() const
    {
        return autoPtr<mesh>(new unstructuredMesh(*this));
    }

    // Destructor
    ~unstructuredMesh();

    //- Disallow default bitwise assignment
    void operator=(const unstructuredMesh&);
};

}

}

#endif
