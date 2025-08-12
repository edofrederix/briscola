#ifndef curvatureScheme_H
#define curvatureScheme_H

#include "dictionary.H"
#include "colocatedFields.H"
#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class twoPhaseModel;

class curvatureScheme
:
    public colocatedScalarField
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the dictionary
    const dictionary& dict_;

    // Reference to the normal scheme
    const normalScheme& normal_;

    // Reference to the void fraction
    const colocatedScalarField& alpha_;

public:

    // Runtime type information
    TypeName("curvatureScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        curvatureScheme,
        dictionary,
        (
            const fvMesh& fvMsh,
            const dictionary& dict,
            const normalScheme& normal,
            const colocatedScalarField& alpha
        ),
        (fvMsh, dict, normal, alpha)
    );

    // Construct from components
    curvatureScheme
    (
        const fvMesh&,
        const dictionary&,
        const normalScheme&,
        const colocatedScalarField&
    );

    // Copy constructor
    curvatureScheme(const curvatureScheme&);

    // Clone normalScheme
    virtual autoPtr<curvatureScheme> clone() const = 0;

    // Selector from components
    static autoPtr<curvatureScheme> New
    (
        const fvMesh&,
        const dictionary&,
        const normalScheme&,
        const colocatedScalarField&
    );

    // Destructor
    virtual ~curvatureScheme();

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the normal scheme
    inline const normalScheme& normal() const { return normal_; }

    // Access to the volume fraction
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Return the curvature field
    inline colocatedScalarField& kappa() { return *this; }

    // Return const access to the curvature field
    inline const colocatedScalarField& kappa() const { return *this; }

    // Compute and return the interpolated kappa field
    tmp<colocatedFaceScalarField> interp() const;

    // Correct the curvature
    virtual void correct() = 0;
};

}

}

}

#endif
