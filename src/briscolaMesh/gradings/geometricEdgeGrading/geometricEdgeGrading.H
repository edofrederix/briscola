#ifndef geometricEdgeGrading_H
#define geometricEdgeGrading_H

#include "grading.H"
#include "tensor.H"
#include "edgeScalar.H"

namespace Foam
{

namespace briscola
{

namespace gradings
{

class geometricEdgeGrading
:
    public grading
{
    // Grading factors
    edgeScalar f_;

    // Disallow default bitwise assignment
    void operator=(const geometricEdgeGrading&);

public:

    // Runtime type information
    TypeName("geometricEdge");

    // Construct from components
    geometricEdgeGrading(const brick& b);

    // Copy constructor
    geometricEdgeGrading(const geometricEdgeGrading&);

    // Clone grading
    virtual autoPtr<grading> clone() const
    {
        return autoPtr<grading>(new geometricEdgeGrading(*this));
    }

    // Linear normalized mapping function
    virtual tmp<scalarField> operator()
    (
        const scalarField& f,
        const label e
    ) const;

    // Transform the grading given a transformation tensor
    virtual void transform(const labelTensor T)
    {
        f_ = pTransform<edgeScalar>(tensor(T), f_);

        for (label cmpt = 0; cmpt < 12; cmpt++)
            if (f_[cmpt] < 0)
                f_[cmpt] = -1.0/f_[cmpt];
    }
};

}

}

}

#endif
