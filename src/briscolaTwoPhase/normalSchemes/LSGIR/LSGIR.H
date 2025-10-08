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
protected:

    // Beta power coefficient (defaults to 1.5)
    scalar beta_;

public:

    // Runtime type information
    TypeName("LSGIR");

    // Construct from components
    LSGIR(const fvMesh&, const dictionary&, const colocatedScalarField&);

    // Copy constructor
    LSGIR(const LSGIR&);

    // Destructor
    virtual ~LSGIR();

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const
    {
        return autoPtr<normalScheme>(new LSGIR(*this));
    }

    // Correct the normal
    virtual void correct();
};

}

}

}

#endif
