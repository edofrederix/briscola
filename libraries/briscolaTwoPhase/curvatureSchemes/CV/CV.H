#ifndef CV_H
#define CV_H

#include "curvatureScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class CV
:
    public curvatureScheme
{

public:

    // Runtime type information
    TypeName("CV");

    // Construct from components
    CV
    (
        const fvMesh&,
        const dictionary&,
        const normalScheme&,
        const colocatedScalarField&
    );

    // Copy constructor
    CV(const CV&);

    // Destructor
    virtual ~CV();

    // Clone curvatureScheme
    virtual autoPtr<curvatureScheme> clone() const
    {
        return autoPtr<curvatureScheme>(new CV(*this));
    }

    // Correct the curvature
    virtual void correct();
};

}

}

}

#endif
