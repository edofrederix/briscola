#ifndef noNormalScheme_H
#define noNormalScheme_H

#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class noNormalScheme
:
    public normalScheme
{
public:

    // Runtime type information
    TypeName("none");

    // Construct from components
    noNormalScheme
    (
        const fvMesh&,
        const dictionary&,
        const colocatedScalarField&
    );

    // Construct with default parameters
    noNormalScheme(const fvMesh&, const colocatedScalarField&);

    // Copy constructor
    noNormalScheme(const noNormalScheme&);

    // Destructor
    virtual ~noNormalScheme();

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const
    {
        return autoPtr<normalScheme>(new noNormalScheme(*this));
    }

    // Correct the normal
    virtual void correct();
};

}

}

}

#endif
