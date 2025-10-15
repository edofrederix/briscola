#ifndef rectilinearMesh_H
#define rectilinearMesh_H

#include "structuredMesh.H"
#include "PartialList.H"
#include "FastPtrList.H"

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
        FastPtrList<scalarList> localCellSizesData_;

        // Local rectilinear points
        FastPtrList<scalarList> localPointsData_;

        // Global rectilinear cell sizes
        FastPtrList<scalarList> globalCellSizesData_;

        // Global rectilinear points
        FastPtrList<scalarList> globalPointsData_;

    // Cell sizes and points without ghosts as partial lists

        // Local rectilinear cell sizes
        FastPtrList<PartialList<scalar>> localCellSizes_;

        // Local rectilinear points
        FastPtrList<PartialList<scalar>> localPoints_;

        // Global rectilinear cell sizes
        FastPtrList<PartialList<scalar>> globalCellSizes_;

        // Global rectilinear points
        FastPtrList<PartialList<scalar>> globalPoints_;

    // Global cell starts for each processor part
    List<labelVector> globalStarts_;

    // Check if global mesh directions are uniform
    labelVector globalUniform_;

    // Set rectilinear metrics
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
    inline const FastPtrList<PartialList<scalar>>& localCellSizes() const
    {
        return localCellSizes_;
    }

    // Local points
    inline const FastPtrList<PartialList<scalar>>& localPoints() const
    {
        return localPoints_;
    }

    // Global cell sizes
    inline const FastPtrList<PartialList<scalar>>& globalCellSizes() const
    {
        return globalCellSizes_;
    }

    // Global points
    inline const FastPtrList<PartialList<scalar>>& globalPoints() const
    {
        return globalPoints_;
    }

    // Return the list of global cell starts for each processor part
    inline const List<labelVector>& globalStarts() const
    {
        return globalStarts_;
    }

    // Return if global mesh is uniform per direction
    inline const labelVector globalUniform() const
    {
        return globalUniform_;
    }

    // Find cells efficiently using the local point data
    virtual labelVector findCell(const vector& p, const label l = 0) const;

    //- Disallow default bitwise assignment
    void operator=(const rectilinearMesh&);
};

}

}

#endif
