#ifndef twoPass_H
#define twoPass_H

#include "CCL.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class twoPass
:
    public CCL
{
    // Labeling stride across processors for parallel tagging
    static const label procStride;

    // Tag union lookup table
    HashTable<label,label> unionTable_;

    // Find root tag
    int find(label);

    // Unite tags
    void uniteTags(label,label);

    // Parallel reduction
    void parallelReduce();

public:

    // Runtime type information
    TypeName("twoPass");

    // Construct from components
    twoPass
    (
        const fvMesh&,
        const dictionary&,
        const colocatedScalarField&
    );

    // Copy constructor
    twoPass(const twoPass&);

    // Destructor
    virtual ~twoPass();

    // Clone CCL
    virtual autoPtr<CCL> clone() const
    {
        return autoPtr<CCL>(new twoPass(*this));
    }

    // Tag the binary field
    virtual void tag();


};

}

}

}

#endif
