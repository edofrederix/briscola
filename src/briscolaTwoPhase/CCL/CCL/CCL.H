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
    colocatedScalarField& alpha_;

    // Number of connected components,
    // reset after each tag() call
    label n_;

    // Minimum size of components (number of cells)
    label minSize_;

    // Total volume of removed alpha
    scalar removedVol_;

    // Tag the binary field
    virtual void tag() = 0;

    // Remove small droplets or bubbles
    void removeSmallComponents();

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
            colocatedScalarField& alpha
        ),
        (fvMsh, dict, alpha)
    );

    // Construct from components
    CCL
    (
        const fvMesh&,
        const dictionary&,
        colocatedScalarField&
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
        colocatedScalarField&
    );

    // Destructor
    virtual ~CCL();

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the volume fraction
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Access the number of connected components
    inline label n() { return n_; }

    // Access to total removed small component volume
    inline scalar removedVol() const { return removedVol_; }

    // Correct the CCL field
    void correct();

    // Assign binary field values
    void binaryTag();

};

}

}

}

#endif
