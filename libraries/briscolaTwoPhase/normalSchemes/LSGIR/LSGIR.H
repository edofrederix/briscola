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

    // Boundary type stencil (true for parallel, periodic and empty boundaries)
    bool boundaryType_[3][3][3];

private:

    // Construct boundaryType
    void createBoundaryType();

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
