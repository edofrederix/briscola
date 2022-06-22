#ifndef symmStencil_H
#define symmStencil_H

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

#include "diagStencil.H"

// The symmStencil class is a compact symmetric scalar stencil, associated with
// the central coefficient and coefficients of the left, bottom and aft neighbor
// of a hex cell. It is defined in the Foam namespace and not in the
// Foam::briscola namespace, just like all other primitives.

namespace Foam
{

class symmStencil;

void writeEntry(Ostream& os, const symmStencil& value);

Istream& operator>>(Istream&, symmStencil&);
Ostream& operator<<(Ostream&, const symmStencil&);

class symmStencil
{

public:

    // The components of this vector space
    scalar v_[4];

    // Component type
    typedef scalar cmptType;

    // Dimensionality of space
    static const direction dim = 3;

    // Number of components in this vector space
    static const direction nComponents = 4;

    //- Rank of symmStencil is 1
    static const direction rank = 1;

    // Component labeling enumeration (left, bottom and aft)
    enum components { C, L, B, A };

    // Static Data Members

    static const char* const typeName;
    static const char* const componentNames[];
    static const symmStencil zero;
    static const symmStencil one;
    static const symmStencil max;
    static const symmStencil min;
    static const symmStencil rootMax;
    static const symmStencil rootMin;

    // Construct null
    inline symmStencil();

    // Construct initialized to zero
    inline symmStencil(const Foam::zero);

    // Construct from four scalars
    inline symmStencil
    (
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&
    );

    // Construct as uniform symmStencil
    inline symmStencil(const scalar&);

    // Construct from Istream
    inline symmStencil(Istream&);

    // Construct as copy of a symmStencil
    inline symmStencil(const symmStencil&);

    // Construct as a copy of a diagStencil
    inline symmStencil(const diagStencil&);

    // Return the number of elements in the symmStencil
    inline static direction size();

    inline const scalar& component(const direction) const;
    inline scalar& component(const direction);

    inline void component(scalar&, const direction) const;
    inline void replace(const direction, const scalar&);

    // Return a symmStencil with all elements being s
    inline static symmStencil uniform(const scalar& s);

    // Access to the center, left, bottom and aft coefficients

    inline const scalar& center() const;
    inline const scalar& left() const;
    inline const scalar& bottom() const;
    inline const scalar& aft() const;

    inline scalar& center();
    inline scalar& left();
    inline scalar& bottom();
    inline scalar& aft();

    // Member Operators

    inline const scalar& operator[](const direction) const;
    inline scalar& operator[](const direction);

    inline void operator+=(const symmStencil&);
    inline void operator+=(const diagStencil&);

    inline void operator-=(const symmStencil&);
    inline void operator-=(const diagStencil&);

    inline void operator=(const scalar);
    inline void operator=(const symmStencil&);
    inline void operator=(const diagStencil&);

    inline void operator=(const Foam::zero);
    inline void operator*=(const scalar);
    inline void operator/=(const scalar);

    // IOstream Operators

    friend Istream& operator>>(Istream&, symmStencil&);
    friend Ostream& operator<<(Ostream&, const symmStencil&);
};

// Return a string representation of a symmStencil
word name(const symmStencil&);

template<>
inline bool contiguous<symmStencil>() {return true;}

}

#include "symmStencilI.H"

#endif
