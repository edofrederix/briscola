#ifndef normalScheme_H
#define normalScheme_H

#include "dictionary.H"
#include "autoPtr.H"
#include "colocatedFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class twoPhaseModel;

class normalScheme
:
    public colocatedVectorField
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the IO dictionary
    const dictionary& dict_;

    // Reference to the volume fraction field
    const colocatedScalarField& alpha_;

public:

    // Runtime type information
    TypeName("normalScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        normalScheme,
        dictionary,
        (
            const fvMesh& fvMsh,
            const dictionary& dict,
            const colocatedScalarField& alpha
        ),
        (fvMsh, dict, alpha)
    );

    // Construct from components
    normalScheme(const fvMesh&, const dictionary&, const colocatedScalarField&);

    // Copy constructor
    normalScheme(const normalScheme&);

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const = 0;

    // Selector from two-phase model
    static autoPtr<normalScheme> New
    (
        const twoPhaseModel&,
        const dictionary&
    );

    // Selector from components
    static autoPtr<normalScheme> New
    (
        const fvMesh&,
        const dictionary&,
        const colocatedScalarField&
    );

    // Destructor
    virtual ~normalScheme();

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the volume fraction field
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Return access to the normal field
    inline colocatedVectorField& normal() { return *this; }

    // Return const access to the normal field
    inline const colocatedVectorField& normal() const { return *this; }

    // Correct the normal
    virtual void correct() = 0;
};

}

}

}

#endif
