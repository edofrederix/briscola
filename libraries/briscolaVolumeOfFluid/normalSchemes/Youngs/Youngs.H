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
    // Normal threshold (defaults to 1e-8)
    scalar threshold_;

public:

    // Runtime type information
    TypeName("Youngs");

    // Construct from components
    Youngs(const vof&, const dictionary&);

    // Copy constructor
    Youngs(const Youngs&);

    // Destructor
    virtual ~Youngs();

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const
    {
        return autoPtr<normalScheme>(new Youngs(*this));
    }

    // Calculate normal
    virtual tmp<colocatedVectorField> operator()();
};

}

}

}

#endif
