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
    // Tensor consisting of rectilinear base vectors on its rows
    tensor base_;

    // Local rectilinear cell sizes
    PtrList<scalarList> localCellSizes_;

    // Local rectlinear points
    PtrList<scalarList> localPoints_;

    // Global rectilinear cell sizes
    PtrList<scalarList> globalCellSizes_;

    // Global rectilinear points
    PtrList<scalarList> globalPoints_;

    // Set rectiliner metrics
    void setMetrics();

public:

    // Runtime type information
    TypeName("rectilinear");

    // Construct from IOdictionary
    rectilinearMesh(const IOdictionary&);

    // Transfer constructor, destroying argument
    rectilinearMesh(autoPtr<mesh>&);

    // Copy constructor
    rectilinearMesh(const rectilinearMesh& msh);

    // Transfer or copy constructor
    rectilinearMesh(rectilinearMesh& msh, bool reuse);

    // Clone mesh
    virtual autoPtr<mesh> clone() const
    {
        return autoPtr<mesh>(new rectilinearMesh(*this));
    }

    // Destructor
    ~rectilinearMesh();

    // Base tensor
    inline tensor base() const { return base_; }

    // Local cell sizes
    inline const PtrList<scalarList>& localCellSizes() const
    {
        return localCellSizes_;
    }

    // Local points
    inline const PtrList<scalarList>& localPoints() const
    {
        return localPoints_;
    }

    // Global cell sizes
    inline const PtrList<scalarList>& globalCellSizes() const
    {
        return globalCellSizes_;
    }

    // Global points
    inline const PtrList<scalarList>& globalPoints() const
    {
        return globalPoints_;
    }

    // Find cells efficiently using the local point data
    virtual labelVector findCell(const vector& p, const label l = 0) const;

    //- Disallow default bitwise assignment
    void operator=(const rectilinearMesh&);
};

}

}

#endif
