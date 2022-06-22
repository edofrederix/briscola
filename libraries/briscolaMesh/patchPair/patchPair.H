#ifndef patchPair_H
#define patchPair_H

#include "patch.H"

namespace Foam
{

namespace briscola
{

class geometry;

class patchPair
:
    public meshObject<geometry>
{
    // Reference to the first patch
    const patch& patch0_;

    // Reference to the second patch
    const patch& patch1_;

    // Check the consistency of this patch pair
    void checkConsistency() const;

public:

    // Construct from geometry, number and two patches
    patchPair(const geometry&, const label, const patch&, const patch&);

    // Copy constructor
    patchPair(const patchPair&);

    // Clone this patchPair (needed by PtrList)
    inline autoPtr<patchPair> clone() const
    {
        return autoPtr<patchPair>(new patchPair(*this));
    }

    // Destructor
    virtual ~patchPair();

    // Access to the first patch
    inline const patch& patch0() const
    {
        return patch0_;
    }

    // Access to the second patch
    inline const patch& patch1() const
    {
        return patch1_;
    }

    // Access to the parent geometry
    inline const geometry& parentGeometry() const
    {
        return parent_;
    }
};

}

}

#endif
