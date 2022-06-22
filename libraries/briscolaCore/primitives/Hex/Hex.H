#ifndef Hex_H
#define Hex_H

#include "direction.H"
#include "scalar.H"
#include "word.H"
#include "zero.H"
#include "contiguous.H"
#include "pTransform.H"
#include "blockUtils.H"

namespace Foam
{

template<class Type>
class Hex;

template<class Type>
void writeEntry(Ostream&, const Hex<Type>&);

template<class Type>
Istream& operator>>(Istream&, Hex<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const Hex<Type>&);

template<class Type>
class Hex
{

public:

    // The components of this vector space
    Type v_[6];

    // Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    // Dimensionality of space
    static const direction dim = Type::dim;

    // Number of components in this vector space
    static const direction nComponents = Type::nComponents;

    //- Rank of type
    static const direction rank = Type::rank;

    // Component labeling enumeration (left, right, bottom, top, aft and fore)
    enum components { L, R, B, T, A, F };

    // Static Data Members

    static const char* const typeName;
    static const char* const componentNames[];
    static const Hex zero;
    static const Hex one;
    static const Hex max;
    static const Hex min;
    static const Hex rootMax;
    static const Hex rootMin;

    // Construct null
    inline Hex();

    // Construct initialized to zero
    inline Hex(const Foam::zero);

    // Construct from six components
    inline Hex
    (
        const Type&,
        const Type&,
        const Type&,
        const Type&,
        const Type&,
        const Type&
    );

    // Construct from Istream
    inline Hex(Istream&);

    // Construct as copy of a Hex
    inline Hex(const Hex<Type>&);

    // Return the number of elements in the Hex
    inline static direction size()
    {
        return 6;
    }

    inline const Type& component(const direction) const;
    inline Type& component(const direction);

    inline void component(Type&, const direction) const;
    inline void replace(const direction, const Type&);

    // Return a Hex with all elements being s
    inline static Hex uniform(const Type& s);

    // Access to left, right, bottom, top, aft and fore element

    inline const Type& left() const;
    inline const Type& right() const;
    inline const Type& bottom() const;
    inline const Type& top() const;
    inline const Type& aft() const;
    inline const Type& fore() const;

    inline Type& left();
    inline Type& right();
    inline Type& bottom();
    inline Type& top();
    inline Type& aft();
    inline Type& fore();

    // Member Operators

    inline Hex<Type> operator+() const;
    inline Hex<Type> operator-() const;

    inline const Type& operator[](const direction) const;
    inline Type& operator[](const direction);

    inline void operator+=(const Hex<Type>&);
    inline void operator-=(const Hex<Type>&);

    inline void operator=(const Foam::zero);
    inline void operator=(const Type&);

    inline void operator*=(const scalar);
    inline void operator*=(const Hex<scalar>&);

    inline void operator/=(const scalar);
    inline void operator/=(const Hex<scalar>&);

    // IOstream Operators

    template<class Type2>
    friend Istream& operator>>(Istream&, Hex<Type2>&);

    template<class Type2>
    friend Ostream& operator<<(Ostream&, const Hex<Type2>&);
};

// Return a string representation of a Hex
template<class Type>
word name(const Hex<Type>&);

template<class Type>
inline Hex<Type> transform
(
    const labelTensor& T,
    const Hex<Type>& h
)
{
    const labelVector x(T & briscola::unitX);
    const labelVector y(T & briscola::unitY);
    const labelVector z(T & briscola::unitZ);

    Hex<Type> h2;

    h2[briscola::faceNumber(-x)] = h.left();
    h2[briscola::faceNumber( x)] = h.right();

    h2[briscola::faceNumber(-y)] = h.bottom();
    h2[briscola::faceNumber( y)] = h.top();

    h2[briscola::faceNumber(-z)] = h.aft();
    h2[briscola::faceNumber( z)] = h.fore();

    return h2;
}

// Hex partial specialization of pTransform. The Hex type has its components
// associated with a local coordinate system, which is subject to transformation
// too.

template<class Type>
class pTransform<Hex<Type>>
{
    Hex<Type> p_;

public:

    // Construct from transformation tensor and primitive
    inline pTransform(const labelTensor& T, const Hex<Type>& p)
    :
        p_(transform(T,p))
    {}

    // Return transformed primitive
    inline operator const Hex<Type>& () const
    {
        return p_;
    }
};

}

#include "HexI.H"

#ifdef NoRepository
    #include "Hex.C"
#endif

#endif
