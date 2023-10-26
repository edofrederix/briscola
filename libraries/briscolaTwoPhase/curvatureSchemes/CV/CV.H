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
    CV(const vof&, const dictionary&);

    // Copy constructor
    CV(const CV&);

    // Destructor
    virtual ~CV();

    // Clone curvatureScheme
    virtual autoPtr<curvatureScheme> clone() const
    {
        return autoPtr<curvatureScheme>(new CV(*this));
    }

    // Calculate curvature
    virtual tmp<colocatedScalarField> operator()(const colocatedVectorField&);
};

}

}

}

#endif
