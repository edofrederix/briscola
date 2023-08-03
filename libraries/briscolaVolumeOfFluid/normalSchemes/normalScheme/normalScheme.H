#ifndef normalScheme_H
#define normalScheme_H

#include "dictionary.H"
#include "colocatedFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class vof;

class normalScheme
{
protected:

    // Reference to the Vof object
    const vof& vf_;

    // Reference to the dictionary
    const dictionary& dict_;

public:

    // Runtime type information
    TypeName("normalScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        normalScheme,
        dictionary,
        (const vof& vf, const dictionary& dict),
        (vf, dict)
    );

    // Construct from components
    normalScheme(const vof&, const dictionary&);

    // Copy constructor
    normalScheme(const normalScheme&);

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const = 0;

    // Selector
    static autoPtr<normalScheme> New(const vof&, const dictionary&);

    // Destructor
    virtual ~normalScheme();

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Compute the normal
    virtual tmp<colocatedVectorField> operator()() = 0;
};

}

}

}

#endif
