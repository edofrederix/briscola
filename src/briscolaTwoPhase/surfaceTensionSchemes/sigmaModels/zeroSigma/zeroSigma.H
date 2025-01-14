#ifndef zeroSigma_H
#define zeroSigma_H

#include "surfaceTensionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class zeroSigma
:
    public surfaceTensionScheme
{
public:

    // Runtime type information
    TypeName("zeroSigma");

    // Construct from components
    zeroSigma
    (
        const twoPhaseModel&,
        const dictionary&
    );

    // Copy constructor
    zeroSigma(const zeroSigma&);

    // Destructor
    virtual ~zeroSigma();

    // Correct sigma
    virtual void correct();
};

}

}

}

#endif
