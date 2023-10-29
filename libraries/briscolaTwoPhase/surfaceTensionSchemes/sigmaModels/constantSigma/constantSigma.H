#ifndef constantSigma_H
#define constantSigma_H

#include "surfaceTensionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class constantSigma
:
    public surfaceTensionScheme
{
public:

    // Runtime type information
    TypeName("constantSigma");

    // Construct from components
    constantSigma
    (
        const fvMesh&,
        const dictionary&,
        const normalScheme&,
        const colocatedScalarField&
    );

    // Copy constructor
    constantSigma(const constantSigma&);

    // Destructor
    virtual ~constantSigma();

    // Correct sigma
    virtual void correct();
};

}

}

}

#endif
