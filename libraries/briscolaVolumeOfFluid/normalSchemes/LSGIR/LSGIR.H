#ifndef LSGIR_H
#define LSGIR_H

#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class LSGIR
:
    public normalScheme
{
    // Normal treshold (defaults to 1e-8)
    scalar threshold_;

public:

    // Runtime type information
    TypeName("LSGIR");

    // Construct from components
    LSGIR(const vof&, const dictionary&);

    // Copy constructor
    LSGIR(const LSGIR&);

    // Destructor
    virtual ~LSGIR();

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const
    {
        return autoPtr<normalScheme>(new LSGIR(*this));
    }

    // Calculate normal
    virtual tmp<colocatedVectorField> operator()();
};

}

}

}

#endif
