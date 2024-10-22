#ifndef immersedBoundaryShape_H
#define immersedBoundaryShape_H

#include "vector.H"
#include "dictionary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class immersedBoundaryShape
{
protected:

    // Inverted bool
    bool inverted_;

public:
    // Runtime type information
    TypeName("immersedBoundaryShape");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        immersedBoundaryShape,
        dictionary,
        (
            const dictionary& dict,
            bool inverted
        ),
        (dict, inverted)
    );

    // Base constructor
    immersedBoundaryShape(const dictionary&, bool);

    // Selector from dictionary
    static autoPtr<immersedBoundaryShape> New(const dictionary&, bool inverted = false);

    // Destructor
    virtual ~immersedBoundaryShape();

    // Check if point is inside shape
    virtual bool isInside(vector) const = 0;

    // Returns dimensionless wall distance given two points The first point is a
    // fluid cell and the second point is an IB cell
    virtual scalar wallDistance(vector, vector) const = 0;

    // Returns the distance from a ghost cell to the immersed boundary in the
    // wall-normal direction
    virtual scalar wallNormalDistance(vector) const = 0;

    // Returns the coordinates of the mirror point of a ghost cell in the
    // wall-normal direction of a shape
    virtual vector mirrorPoint(vector) const = 0;

    // Inverted shape
    virtual bool inverted() const { return inverted_; }
};

}

}

}

#endif
