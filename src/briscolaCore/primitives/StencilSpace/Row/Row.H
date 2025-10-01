#ifndef Row_H
#define Row_H

#include "FixedList.H"

// Contiguous container of row data, putting everything in a single scalar list.
// This is useful for parallel data transfer in one call. Based on a fixed list
// of size equal to the number of stencil elements plus the number of source
// type components.

namespace Foam
{

template<class SType, class Type>
class Row
:
    public FixedList
    <
        scalar,
        SType::nCsComponents
      + pTraits<Type>::nComponents
    >
{
public:

    typedef FixedList
    <
        scalar,
        SType::nCsComponents+pTraits<Type>::nComponents
    > listType;

    // Null constructor
    inline Row();

    // Construct from components
    inline Row(const SType& s, const Type& b);

    // Return the stencil
    inline SType stencil() const;

    // Return the source
    inline Type source() const;
};

}

#include "RowI.H"

#endif