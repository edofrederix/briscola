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
    // Disallow default bitwise assignment
    void operator=(const Youngs&);

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
