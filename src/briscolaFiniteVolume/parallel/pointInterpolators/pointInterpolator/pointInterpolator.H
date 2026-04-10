#ifndef pointInterpolator_H
#define pointInterpolator_H

#include "fvMesh.H"
#include "meshField.H"
#include "SubList.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class pointInterpolator
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Is the provided list of points global? If it is global, it means all
    // processors have a copy of all points. If it is local, each processor has
    // its own set of points.
    const Switch global_;

    // Size of this processor's part, only relevant when the provided list of
    // points is local
    label size_;

    // Mesh level
    label l_;

    // Mesh direction
    label d_;

    // List of points to interpolate to. The first size_ points are ours.
    // Additional points in the list are others'.
    DynamicList<vector> points_;

    // Requester processor of each point
    DynamicList<label> requesters_;

    // Provider processor of each point
    DynamicList<label> providers_;

    // Cell indices belonging to the points
    DynamicList<labelVector> indices_;

    // Cell coordinates of the points
    DynamicList<vector> cellCoordinates_;

    // Gather and scatter data
    template<class Type>
    List<Type> combine(List<Type>&) const;

private:

    // Disallow default bitwise assignment
    void operator=(const pointInterpolator&);

public:

    // Runtime type information
    TypeName("pointInterpolator");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        pointInterpolator,
        dictionary,
        (
            const fvMesh& fvMsh,
            const vectorList& points,
            const bool global,
            const label l,
            const label d
        ),
        (fvMsh, points, global, l, d)
    );

    // Construct from mesh and points
    pointInterpolator
    (
        const fvMesh& fvMsh,
        const vectorList& points,
        const bool global = true,
        const label l = 0,
        const label d = 0
    );

    // Copy constructor
    pointInterpolator(const pointInterpolator<MeshType>&);

    // Selector from a list of points that can be global or local
    static autoPtr<pointInterpolator<MeshType>> New
    (
        const fvMesh& fvMsh,
        const vectorList& points,
        const word type,
        const bool global = true,
        const label l = 0,
        const label d = 0
    );

    // Clone the interpolator
    virtual autoPtr<pointInterpolator<MeshType>> clone() const = 0;

    // Destructor
    virtual ~pointInterpolator();

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return the mesh level
    inline label l() const { return l_; }

    // Return the mesh direction
    inline label d() const { return d_; }

    // Number of points
    inline label size() const { return size_; }

    // Return points as sub-list
    inline SubList<vector> points() const
    {
        return SubList<vector>(points_, size_);
    }

    // Return list of missing points
    virtual vectorList missingPoints() const;

    // Interpolators

    #define INTERPFUNC(TYPE) \
    virtual List<TYPE> operator() (const meshField<TYPE,MeshType>&) const = 0;

    INTERPFUNC(scalar);
    INTERPFUNC(vector);
    INTERPFUNC(tensor);
    INTERPFUNC(sphericalTensor);
    INTERPFUNC(symmTensor);
    INTERPFUNC(diagTensor);
    INTERPFUNC(faceScalar);
    INTERPFUNC(faceVector);

    #undef INTERPFUNC
};

}

}

}

#include "pointInterpolatorTemplates.C"

#endif
