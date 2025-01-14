#ifndef MYC_H
#define MYC_H

#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class MYC
:
    public normalScheme
{
public:

    // Runtime type information
    TypeName("MYC");

    // Construct from components
    MYC(const fvMesh&, const dictionary&, const colocatedScalarField&);

    // Copy constructor
    MYC(const MYC&);

    // Destructor
    virtual ~MYC();

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const
    {
        return autoPtr<normalScheme>(new MYC(*this));
    }

    // Correct the normal
    virtual void correct();
};

}

}

}

#endif
