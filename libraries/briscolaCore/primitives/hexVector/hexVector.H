#ifndef hexVector_H
#define hexVector_H

#include "Hex.H"
#include "vector.H"

namespace Foam
{

typedef Hex<vector> hexVector;

template<>
inline bool contiguous<hexVector>() {return true;}

template<>
class pTraits<hexVector>
{
    hexVector p_;

public:

    // Component type
    typedef scalar cmptType;

    // Equivalent type of labels used for valid component indexing
    typedef label labelType;

    // Dimensionality of space
    static const direction dim = 3;

    // Rank of hexVector is 1
    static const direction rank = 1;

    // Number of components in hexVector is 6
    static const direction nComponents = 6;

    // Static Data Members

    static const char* const typeName;
    static const char* const componentNames[];
    static const hexVector zero;
    static const hexVector one;
    static const hexVector max;
    static const hexVector min;
    static const hexVector rootMax;
    static const hexVector rootMin;

    // Construct from primitive
    explicit pTraits(const hexVector&);

    // Construct from Istream
    pTraits(Istream&);

    // Access to the hexVector value
    operator hexVector() const
    {
        return p_;
    }
};

}

#endif
