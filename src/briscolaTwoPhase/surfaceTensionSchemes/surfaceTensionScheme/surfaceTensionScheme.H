#ifndef surfaceTensionScheme_H
#define surfaceTensionScheme_H

#include "dictionary.H"
#include "colocatedFields.H"
#include "normalScheme.H"
#include "curvatureScheme.H"
#include "interpolationScheme.H"

// The surface tension scheme computes the surface tension potential, defined as
// a flux in the form phi = A_f*(kappa*sigma)_f*grad_f(alpha).

namespace Foam
{

namespace briscola
{

namespace fv
{

class twoPhaseModel;

class surfaceTensionScheme
:
    public colocatedScalarFaceField
{
protected:

    // Reference to the two-phase model
    const twoPhaseModel& tpm_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the dictionary
    const dictionary& dict_;

    // Reference to the normal scheme
    const normalScheme& normal_;

    // Reference to the void fraction
    const colocatedScalarField& alpha_;

    // Surface tension field
    colocatedScalarField sigma_;

    // Curvature scheme
    autoPtr<curvatureScheme> curvatureSchemePtr_;

public:

    // Runtime type information
    TypeName("surfaceTensionScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        surfaceTensionScheme,
        dictionary,
        (
            const twoPhaseModel& tpm,
            const dictionary& dict
        ),
        (tpm, dict)
    );

    // Construct from two-phase model
    surfaceTensionScheme
    (
        const twoPhaseModel&,
        const dictionary&
    );

    // Copy constructor
    surfaceTensionScheme(const surfaceTensionScheme&);

    // Clone surfaceTensionScheme
    virtual autoPtr<surfaceTensionScheme> clone() const = 0;

    // Selector from two-phase model
    static autoPtr<surfaceTensionScheme> New
    (
        const twoPhaseModel&,
        const dictionary&
    );

    // Destructor
    virtual ~surfaceTensionScheme();

    // Access to the two-phase model
    inline const twoPhaseModel& tpm() const { return tpm_; }

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the normal scheme
    inline const normalScheme& normal() const { return normal_; }

    // Access to the volume fraction field
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Access to the curvature scheme
    inline const curvatureScheme& kappa() const
    {
        return curvatureSchemePtr_();
    }

    // Access to the sigma field
    inline const colocatedScalarField& sigma() const { return sigma_; }

    // Correct the surface tension model
    virtual void correct();
};

}

}

}

#endif
