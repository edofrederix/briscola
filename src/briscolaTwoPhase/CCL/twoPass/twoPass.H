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
private:

    // Tag union lookup table
    HashTable<label,label> unionTable_;

    // Find root tag in union table
    int find(label);

    // Find root tag in table
    int find(label, HashTable<label,label>&);

    // Unite tags
    void uniteTags(label,label);

    // Parallel reduction
    void parallelReduce();

protected:

    // Tag the binary field
    virtual void tag();

public:

    // Runtime type information
    TypeName("twoPass");

    // Construct from components
    twoPass
    (
        const fvMesh&,
        const dictionary&,
        colocatedScalarField&
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

};

}

}

}

#endif
