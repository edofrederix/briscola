#ifndef periodicPatch_H
#define periodicPatch_H

#include "patch.H"

namespace Foam
{

namespace briscola
{

class periodicPatch
:
    public patch
{
private:

    // Name of the neighbor
    word neighbor_;

public:

    // Runtime type information
    TypeName("periodic");

    // Construct from dictionary
    periodicPatch(const geometry&, const label, const word, const dictionary&);

    // Copy constructor
    periodicPatch(const periodicPatch&);

    // Clone this periodic patch
    inline autoPtr<periodicPatch> clone() const
    {
        return autoPtr<periodicPatch>(new periodicPatch(*this));
    }

    // Destructor
    virtual ~periodicPatch();

    // Name of the neighbor
    inline word neighbor() const
    {
        return neighbor_;
    }
};

}

}

#endif
