#ifndef curvatureScheme_H
#define curvatureScheme_H

#include "dictionary.H"
#include "colocatedFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class vof;

class curvatureScheme
{
protected:

    // Reference to the Vof object
    const vof& vf_;

    // Reference to the dictionary
    const dictionary& dict_;

public:

    // Runtime type information
    TypeName("curvatureScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        curvatureScheme,
        dictionary,
        (const vof& vf, const dictionary& dict),
        (vf, dict)
    );

    // Construct from components
    curvatureScheme(const vof&, const dictionary&);

    // Copy constructor
    curvatureScheme(const curvatureScheme&);

    // Clone normalScheme
    virtual autoPtr<curvatureScheme> clone() const = 0;

    // Selector
    static autoPtr<curvatureScheme> New(const vof&, const dictionary&);

    // Destructor
    virtual ~curvatureScheme();

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Compute the curvature
    virtual tmp<colocatedScalarField> operator()(const colocatedVectorField&) = 0;
};

}

}

}

#endif
