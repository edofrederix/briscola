#ifndef cachedRefCount_H
#define cachedRefCount_H

#include "refCount.H"

namespace Foam
{

class cachedRefCount
:
    public refCount
{

    // Is the object that contains this counter cached?
    bool cached_;

public:

    // Construct null with cached bool set to false
    cachedRefCount()
    :
        cached_(false)
    {}

    // Return the reference count subtracted by one if the object is cached,
    // to allow for an additional reference from the cache registry
    int count() const
    {
        return refCount::count() - label(cached_);
    }

    // Access to the cached bool
    bool& cached() { return cached_; }
};

}

#endif
