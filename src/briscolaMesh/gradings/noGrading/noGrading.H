#ifndef noGrading_H
#define noGrading_H

#include "grading.H"

namespace Foam
{

namespace briscola
{

namespace gradings
{

class noGrading
:
    public grading
{
    // Disallow default bitwise assignment
    void operator=(const noGrading&);

public:

    // Runtime type information
    TypeName("none");

    // Construct from components
    noGrading(const brick& b);

    // Copy constructor
    noGrading(const noGrading&);

    // Clone grading
    virtual autoPtr<grading> clone() const
    {
        return autoPtr<grading>(new noGrading(*this));
    }

    // Linear normalized mapping function
    virtual tmp<scalarField> operator()
    (
        const scalarField& f,
        const label e
    ) const
    {
        return tmp<scalarField>(new scalarField(f));
    }

    // Transform the grading given a transformation tensor
    virtual void transform(const labelTensor)
    {}
};

}

}

}

#endif
