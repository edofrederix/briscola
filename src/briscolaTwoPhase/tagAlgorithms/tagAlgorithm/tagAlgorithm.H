#ifndef tagAlgorithm_H
#define tagAlgorithm_H

#include "dictionary.H"
#include "colocatedFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class tagAlgorithm
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

    // Number of particles,
    // reset after each tag() call
    label n_;

    // Minimum size of particles (number of cells)
    label minSize_;

    // Total volume of removed alpha
    scalar removedVol_;

    // Tag the binary field
    virtual void tag() = 0;

    // Remove small particles
    void removeSmallParticles();

public:

    // Runtime type information
    TypeName("tagAlgorithm");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        tagAlgorithm,
        dictionary,
        (
            const fvMesh& fvMsh,
            const dictionary& dict,
            colocatedScalarField& alpha
        ),
        (fvMsh, dict, alpha)
    );

    // Construct from components
    tagAlgorithm
    (
        const fvMesh&,
        const dictionary&,
        colocatedScalarField&
    );

    // Copy constructor
    tagAlgorithm(const tagAlgorithm&);

    // Clone tag algorithm
    virtual autoPtr<tagAlgorithm> clone() const = 0;

    // Selector from components
    static autoPtr<tagAlgorithm> New
    (
        const fvMesh&,
        const dictionary&,
        colocatedScalarField&
    );

    // Destructor
    virtual ~tagAlgorithm();

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the volume fraction
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Access the number of particles
    inline label n() { return n_; }

    // Access to total removed small particles volume
    inline scalar removedVol() const { return removedVol_; }

    // Correct the tag field
    void correct();

    // Assign binary field values
    void binaryTag();

};

}

}

}

#endif
