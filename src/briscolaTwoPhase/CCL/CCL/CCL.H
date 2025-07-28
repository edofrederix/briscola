#ifndef CCL_H
#define CCL_H

#include "dictionary.H"
#include "colocatedFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class CCL
:
    public colocatedLabelField
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the dictionary
    const dictionary& dict_;

    // Reference to the void fraction
    const colocatedScalarField& alpha_;

public:

    // Runtime type information
    TypeName("CCL");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        CCL,
        dictionary,
        (
            const fvMesh& fvMsh,
            const dictionary& dict,
            const colocatedScalarField& alpha
        ),
        (fvMsh, dict, alpha)
    );

    // Construct from components
    CCL
    (
        const fvMesh&,
        const dictionary&,
        const colocatedScalarField&
    );

    // Copy constructor
    CCL(const CCL&);

    // Clone normalScheme
    virtual autoPtr<CCL> clone() const = 0;

    // Selector from components
    static autoPtr<CCL> New
    (
        const fvMesh&,
        const dictionary&,
        const colocatedScalarField&
    );

    // Destructor
    virtual ~CCL();

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the volume fraction
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Assign binary field values
    void binaryTag();

    // Tag the binary field
    virtual void tag() = 0;
};

}

}

}

#endif
