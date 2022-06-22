#ifndef hexScalar_H
#define hexScalar_H

#include "Hex.H"
#include "scalar.H"

namespace Foam
{

typedef Hex<scalar> hexScalar;

template<>
inline bool contiguous<hexScalar>() {return true;}

template<>
class pTraits<hexScalar>
{
    hexScalar p_;

public:

    // Component type
    typedef scalar cmptType;

    // Equivalent type of labels used for valid component indexing
    typedef label labelType;

    // Dimensionality of space
    static const direction dim = 3;

    // Rank of hexScalar is 0
    static const direction rank = 1;

    // Number of components in hexScalar is 6
    static const direction nComponents = 6;

    // Static Data Members

    static const char* const typeName;
    static const char* const componentNames[];
    static const hexScalar zero;
    static const hexScalar one;
    static const hexScalar max;
    static const hexScalar min;
    static const hexScalar rootMax;
    static const hexScalar rootMin;

    // Construct from primitive
    explicit pTraits(const hexScalar&);

    // Construct from Istream
    pTraits(Istream&);

    // Access to the hexScalar value
    operator hexScalar() const
    {
        return p_;
    }
};

}

#endif
