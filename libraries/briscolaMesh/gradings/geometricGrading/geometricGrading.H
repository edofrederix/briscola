#ifndef geometricGrading_H
#define geometricGrading_H

#include "grading.H"

namespace Foam
{

namespace briscola
{

namespace gradings
{

class geometricGrading
:
    public grading
{
    // Grading factors
    vector f_;

    // Copy constructor
    geometricGrading(const geometricGrading&);

    // Disallow default bitwise assignment
    void operator=(const geometricGrading&);

public:

    // Runtime type information
    TypeName("geometric");

    // Construct from components
    geometricGrading(const brick& b);

    // Linear normalized mapping function
    virtual tmp<scalarField> operator()
    (
        const scalarField& f,
        const label e
    ) const;

    // Transform the grading given a transformation tensor
    virtual void transform(const labelTensor T)
    {
        f_ = tensor(T) & f_;

        for (label cmpt = 0; cmpt < 3; cmpt++)
            if (f_[cmpt] < 0)
                f_[cmpt] = -1.0/f_[cmpt];
    }
};

}

}

}

#endif
