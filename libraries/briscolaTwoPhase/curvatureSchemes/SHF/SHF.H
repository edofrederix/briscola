#ifndef SHF_H
#define SHF_H

#include "curvatureScheme.H"
#include "cellDataExchange.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class SHF
:
    public curvatureScheme
{

    // Cell data exchange object
    autoPtr<cellDataExchange<colocated>> exchangePtr_;

    // Exchange data sizes per face
    labelList exchangeSizes_;

    // Initialize communication
    void setCommunication();

public:

    // Runtime type information
    TypeName("SHF");

    // Construct from components
    SHF
    (
        const fvMesh&,
        const dictionary&,
        const normalScheme&,
        const colocatedScalarField&
    );

    // Copy constructor
    SHF(const SHF&);

    // Destructor
    virtual ~SHF();

    // Clone curvatureScheme
    virtual autoPtr<curvatureScheme> clone() const
    {
        return autoPtr<curvatureScheme>(new SHF(*this));
    }

    // Correct the curvature
    virtual void correct();
};

}

}

}

#endif
