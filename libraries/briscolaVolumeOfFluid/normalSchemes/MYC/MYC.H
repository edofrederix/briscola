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
    // Normal treshold (defaults to 1e-8)
    scalar threshold_;

public:

    // Runtime type information
    TypeName("MYC");

    // Construct from components
    MYC(const vof&, const dictionary&);

    // Copy constructor
    MYC(const MYC&);

    // Destructor
    virtual ~MYC();

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const
    {
        return autoPtr<normalScheme>(new MYC(*this));
    }

    // Calculate normal
    virtual tmp<colocatedVectorField> operator()();
};

}

}

}

#endif
