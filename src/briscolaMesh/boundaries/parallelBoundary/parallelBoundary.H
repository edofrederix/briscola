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

    // Rotation of the neighbor in case of a parallel boundary (defaults to the
    // identity tensor)
    labelTensor T_;

    // Master bool. Boundaries that require communication can be master or
    // slave. By default a boundary is master. Can be modified by derived
    // classes.
    bool master_;

    // Communication tag
    label tag_;

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

    // Rotation of the neighbor
    inline const labelTensor& T() const { return T_; }

    // Master bool
    inline bool master() const { return master_; }
    inline bool slave() const { return !master_; }

    // Tag number
    virtual label tag() const
    {
        return tag_;
    }

    // Set the tag number
    virtual void setTag(const label t)
    {
        tag_ = t;
    }

    // Neighbor processor number
    virtual label neighborProcNum() const { return neighborProcNum_; }
};

}

}

#endif
