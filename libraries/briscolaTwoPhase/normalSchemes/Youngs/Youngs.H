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
    // Volume fraction gradient threshold (defaults to 1e-8)
    scalar gradThreshold_;

public:

    // Runtime type information
    TypeName("Youngs");

    // Construct from components
    Youngs(const fvMesh&, const dictionary&, const colocatedScalarField&);

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
