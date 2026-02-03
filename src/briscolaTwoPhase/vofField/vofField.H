#ifndef vofField_H
#define vofField_H

#include "dictionary.H"
#include "colocatedFields.H"
#include "normalScheme.H"
#include "vof.H"
#include "surfaceTensionScheme.H"
#include "tagAlgorithm.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class vofField
:
    public colocatedScalarField
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the dictionary
    const dictionary& dict_;

    // Normal scheme
    autoPtr<normalScheme> normalSchemePtr_;

    // Surface tension scheme
    autoPtr<surfaceTensionScheme> surfaceTensionSchemePtr_;

    // The Vof model
    autoPtr<vof> vfPtr_;

    // Tag algorithm
    autoPtr<tagAlgorithm> tagAlgorithmPtr_;

public:

    // Runtime type information
    TypeName("vofField");

    // Construct from components
    vofField
    (
        const colocatedScalarField&,
        const dictionary&,
        bool initToZero = false
    );

    // Copy constructor
    vofField(const vofField&);

    // Clone vofField
    virtual autoPtr<vofField> clone() const
    {
        return autoPtr<vofField>(new vofField(*this));
    }

    // Destructor
    virtual ~vofField();

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Return access to the normal scheme
    inline normalScheme& normal() { return normalSchemePtr_(); }

    // Return access to the surface tension sheme
    inline surfaceTensionScheme& surfaceTension() { return surfaceTensionSchemePtr_(); }

    // Return access to the Vof model
    inline vof& vf() { return vfPtr_(); }

    // Return access to the tag algorithm
    inline tagAlgorithm& tag() { return tagAlgorithmPtr_(); }

};

}

}

}

#endif
