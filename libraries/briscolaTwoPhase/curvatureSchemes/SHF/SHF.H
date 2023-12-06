#ifndef SHF_H
#define SHF_H

#include "curvatureScheme.H"

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

    // List of the indices outside the processor for
    // each face boundary
    List<List<labelVector>> communicationIndices_;

    // List of the kappa values outside the processor
    // for each face boundary
    List<List<scalar>> commKappa_;

    // Maximum allowed value of kappa
    scalar maxKappa_;

    // Initialize communication indices
    void setCommunicationIndices();

    // Initialize maximum allowed value of kappa
    void setMaxKappa();

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
