#ifndef Youngs_H
#define Youngs_H

#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class Youngs
:
    public normalScheme
{
    // Number of alpha smoothing cycles (defaults to 1)
    label nSmooth_;

public:

    // Runtime type information
    TypeName("Youngs");

    // Construct from components
    Youngs(const fvMesh&, const dictionary&, const colocatedScalarField&);

    // Construct with default parameters
    Youngs(const fvMesh&, const colocatedScalarField&);

    // Copy constructor
    Youngs(const Youngs&);

    // Destructor
    virtual ~Youngs();

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const
    {
        return autoPtr<normalScheme>(new Youngs(*this));
    }

    // Correct the normal
    virtual void correct();
};

}

}

}

#endif
