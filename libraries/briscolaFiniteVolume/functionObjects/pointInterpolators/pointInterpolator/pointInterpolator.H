#ifndef pointInterpolator_H
#define pointInterpolator_H

#include "fvMesh.H"
#include "meshField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class pointInterpolator
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // List of points to interpolate to
    vectorList points_;

    // Cell indices belonging to the points
    List<labelVector> indices_;

    // Cell coordinates of the points
    vectorList cellCoordinates_;

    // Gather and scatter data
    template<class Type>
    void gatherScatter(List<Type>&);

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
            const vectorList& points
        ),
        (fvMsh, points)
    );

    // Construct from mesh and points
    pointInterpolator(const fvMesh& fvMsh, const vectorList& points);

    // Copy constructor
    pointInterpolator(const pointInterpolator&);

    // Selector
    static autoPtr<pointInterpolator> New
    (
        const fvMesh& fvMsh,
        const vectorList& point,
        const word type
    );

    // Destructor
    virtual ~pointInterpolator();

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the points list
    inline const vectorList& points() const { return points_; }

    // Access to the cell indices
    inline const List<labelVector>& indices() const { return indices_; }

    // Access to the cell coordinates
    inline const vectorList& cellCoordinates() const
    {
        return cellCoordinates_;
    }

    // Check if all points were found
    virtual bool allPointsFound() const;

    // Interpolators

    #define INTERPFUNC(TYPE) \
    virtual List<TYPE> operator() (const meshField<TYPE,colocated>&) = 0;

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
