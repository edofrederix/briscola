#ifndef PartialList_H
#define PartialList_H

#include "SubList.H"
#include "autoPtr.H"

namespace Foam
{

// Copy of SubList, but then with an alternative []-operator to allow for
// negative indexing.

template<class T>
class PartialList
:
    public SubList<T>
{
    // Size of underlying data
    label dataSize_;

    // Start index of partial list
    label startIndex_;

    // Check index i is within valid data range (0 ... size-1)
    inline void checkIndex(const label) const;

public:

    // Return a null PartialList
    inline static const PartialList<T>& null();

    // Construct from UList and sub-list size
    inline PartialList
    (
        const UList<T>& list,
        const label subSize
    );

    // Construct from UList, sub-list size and start index
    inline PartialList
    (
        const UList<T>& list,
        const label subSize,
        const label startIndex
    );

    // Copy constructor
    inline PartialList(const PartialList<T>& sl) = default;

    // Clone
    inline autoPtr<PartialList<T>> clone() const;

    // Return element of PartialList
    inline T& operator[](const label);

    // Return const element of PartialList
    inline const T& operator[](const label) const;
};


}

#include "PartialListI.H"

#endif
