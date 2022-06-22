#ifndef stencil_H
#define stencil_H

#include "direction.H"
#include "scalar.H"
#include "vector.H"
#include "diagTensor.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "tensor.H"
#include "word.H"
#include "zero.H"
#include "contiguous.H"

#include "symmStencil.H"
#include "diagStencil.H"
#include "hexScalar.H"

// The stencil class is a compact scalar stencil, associated with the central
// coefficient and the coefficients of the left, right, bottom, top, aft and
// fore neighbor of a hex cell. It is defined in the Foam namespace and not in
// the Foam::briscola namespace, just like all other primitives.

namespace Foam
{

class stencil;

void writeEntry(Ostream& os, const stencil& value);

Istream& operator>>(Istream&, stencil&);
Ostream& operator<<(Ostream&, const stencil&);

class stencil
{

public:

    // The components of this vector space
    scalar v_[7];

    // Component type
    typedef scalar cmptType;

    // Dimensionality of space
    static const direction dim = 3;

    // Number of components in this vector space
    static const direction nComponents = 7;

    //- Rank of stencil is 1
    static const direction rank = 1;

    // Component labeling enumeration (central, left, rigt, bottom, top, aft and
    // fore)
    enum components { C, L, R, B, T, A, F };

    // Static Data Members

    static const char* const typeName;
    static const char* const componentNames[];
    static const stencil zero;
    static const stencil one;
    static const stencil max;
    static const stencil min;
    static const stencil rootMax;
    static const stencil rootMin;

    // Construct null
    inline stencil();

    // Construct initialized to zero
    inline stencil(const Foam::zero);

    // Construct from seven scalars
    inline stencil
    (
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&
    );

    // Construct as uniform stencil
    inline stencil(const scalar&);

    // Construct from Istream
    inline stencil(Istream&);

    // Construct as copy of a stencil
    inline stencil(const stencil&);

    // Construct as a copy of a symmStencil
    inline stencil(const symmStencil&);

    // Construct as a copy of a diagStencil
    inline stencil(const diagStencil&);

    // Construct as a copy of a hexScalar
    inline stencil(const hexScalar&);

    // Return the number of elements in the stencil
    inline static direction size();

    inline const scalar& component(const direction) const;
    inline scalar& component(const direction);

    inline void component(scalar&, const direction) const;
    inline void replace(const direction, const scalar&);

    // Return a stencil with all elements being s
    inline static stencil uniform(const scalar& s);

    // Access to center, left, right, bottom, top, aft and fore coefficients

    inline const scalar& center() const;
    inline const scalar& left() const;
    inline const scalar& right() const;
    inline const scalar& bottom() const;
    inline const scalar& top() const;
    inline const scalar& aft() const;
    inline const scalar& fore() const;

    inline scalar& center();
    inline scalar& left();
    inline scalar& right();
    inline scalar& bottom();
    inline scalar& top();
    inline scalar& aft();
    inline scalar& fore();

    // Member Operators

    inline const scalar& operator[](const direction) const;
    inline scalar& operator[](const direction);

    inline void operator+=(const stencil&);
    inline void operator+=(const symmStencil&);
    inline void operator+=(const diagStencil&);

    inline void operator-=(const stencil&);
    inline void operator-=(const symmStencil&);
    inline void operator-=(const diagStencil&);

    inline void operator=(const scalar);
    inline void operator=(const stencil&);
    inline void operator=(const symmStencil&);
    inline void operator=(const diagStencil&);
    inline void operator=(const hexScalar&);

    inline void operator=(const Foam::zero);
    inline void operator*=(const scalar);
    inline void operator/=(const scalar);

    // IOstream Operators

    friend Istream& operator>>(Istream&, stencil&);
    friend Ostream& operator<<(Ostream&, const stencil&);
};

// Return a string representation of a stencil
word name(const stencil&);

template<>
inline bool contiguous<stencil>() {return true;}

}

#include "stencilI.H"

#endif
