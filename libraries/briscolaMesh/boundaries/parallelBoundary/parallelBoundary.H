#ifndef parallelBoundary_H
#define parallelBoundary_H

#include "boundary.H"

namespace Foam
{

namespace briscola
{

class parallelBoundary
:
    public boundary
{
    // Processor number of the other side
    label neighborProcNum_;

    // Disallow default bitwise assignment
    void operator=(const parallelBoundary&);

public:

    // Type number
    const static label typeNumber;

    // Runtime type information
    TypeName("parallel");

    // Construct from components
    parallelBoundary(const mesh& msh, const dictionary& dict);

    // Copy constructor
    parallelBoundary(const parallelBoundary&);

    // Clone boundary
    virtual autoPtr<boundary> clone() const
    {
        return autoPtr<boundary>(new parallelBoundary(*this));
    }

    // Boundary type number
    virtual label typeNum() const
    {
        return typeNumber;
    }

    // Neighbor processor number
    virtual label neighborProcNum() const { return neighborProcNum_; }
};

}

}

#endif
