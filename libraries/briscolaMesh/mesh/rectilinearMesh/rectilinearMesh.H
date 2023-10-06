#ifndef rectilinearMesh_H
#define rectilinearMesh_H

#include "structuredMesh.H"
#include "PartialList.H"

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

    // Cell sizes and points with ghosts

        // Local rectilinear cell sizes
        PtrList<scalarList> localCellSizesData_;

        // Local rectilinear points
        PtrList<scalarList> localPointsData_;

        // Global rectilinear cell sizes
        PtrList<scalarList> globalCellSizesData_;

        // Global rectilinear points
        PtrList<scalarList> globalPointsData_;

    // Cell sizes and points without ghosts as partial lists

        // Local rectilinear cell sizes
        PtrList<PartialList<scalar>> localCellSizes_;

        // Local rectilinear points
        PtrList<PartialList<scalar>> localPoints_;

        // Global rectilinear cell sizes
        PtrList<PartialList<scalar>> globalCellSizes_;

        // Global rectilinear points
        PtrList<PartialList<scalar>> globalPoints_;

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
    inline const tensor& base() const { return base_; }

    // Cell sizes and points without ghosts as partial lists. Ghost cell
    // indexing (at -1 or N) is allowed.

    // Local cell sizes
    inline const PtrList<PartialList<scalar>>& localCellSizes() const
    {
        return localCellSizes_;
    }

    // Local points
    inline const PtrList<PartialList<scalar>>& localPoints() const
    {
        return localPoints_;
    }

    // Global cell sizes
    inline const PtrList<PartialList<scalar>>& globalCellSizes() const
    {
        return globalCellSizes_;
    }

    // Global points
    inline const PtrList<PartialList<scalar>>& globalPoints() const
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
