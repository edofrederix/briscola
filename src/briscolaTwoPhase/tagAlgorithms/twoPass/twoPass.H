#ifndef twoPass_H
#define twoPass_H

#include "tagAlgorithm.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class twoPass
:
    public tagAlgorithm
{
private:

    // Alias for HashTable needed for forAllConstIter macro
    typedef HashTable<label,label> labelTable;

    // Tag union lookup table
    HashTable<label,label> unionTable_;

    // Find root tag in union table
    label find(const label);

    // Find root tag in table
    label find(const label, HashTable<label,label>&);

    // Unite tags in union table
    void uniteTags(const label, const label);

    // Unite tags in table
    void uniteTags(HashTable<label,label>&, const label, const label);

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

    // Clone tagAlgorithm
    virtual autoPtr<tagAlgorithm> clone() const
    {
        return autoPtr<tagAlgorithm>(new twoPass(*this));
    }

};

}

}

}

#endif
