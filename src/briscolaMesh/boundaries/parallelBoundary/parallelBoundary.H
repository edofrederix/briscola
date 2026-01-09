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

    // Communication tag
    label tag_;

    // Disallow default bitwise assignment
    void operator=(const parallelBoundary&);

public:

    // Runtime type information
    TypeName("parallel");

    // Construct from components
    parallelBoundary(const level& lvl, const dictionary& dict);

    // Copy constructor
    parallelBoundary(const parallelBoundary&);

    // Copy constructor with alternative level reference
    parallelBoundary(const parallelBoundary&, const level&);

    // Clone boundary
    virtual autoPtr<boundary> clone() const
    {
        return autoPtr<boundary>(new parallelBoundary(*this));
    }

    // Clone boundary with alternative level reference
    virtual autoPtr<boundary> clone(const level& lvl) const
    {
        return autoPtr<boundary>(new parallelBoundary(*this, lvl));
    }

    // Rotation of the neighbor
    inline const labelTensor& T() const { return T_; }

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
