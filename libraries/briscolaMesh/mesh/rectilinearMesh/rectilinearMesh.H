#ifndef rectilinearMesh_H
#define rectilinearMesh_H

#include "structuredMesh.H"

namespace Foam
{

namespace briscola
{

class rectilinearMesh
:
    public structuredMesh
{
    // Rectilinear cell sizes
    PtrList<scalarList> cellSizes_;

    // Set rectilinear cell sizes
    void setCellSizes();

public:

    // Runtime type information
    TypeName("rectilinear");

    // Construct from IOdictionary
    rectilinearMesh(const IOdictionary&);

    // Transfer constructor, destroying argument
    rectilinearMesh(autoPtr<mesh>&);

    // Transfer or copy constructor
    rectilinearMesh(rectilinearMesh& msh, bool reuse);

    // Clone mesh
    virtual autoPtr<mesh> clone() const
    {
        return autoPtr<mesh>(new rectilinearMesh(*this));
    }

    // Destructor
    ~rectilinearMesh();

    inline const PtrList<scalarList>& cellSizes() const
    {
        return cellSizes_;
    }

    //- Disallow default bitwise assignment
    void operator=(const rectilinearMesh&);
};

}

}

#endif
