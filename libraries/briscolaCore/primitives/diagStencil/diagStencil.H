#ifndef diagStencil_H
#define diagStencil_H

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

// The diagStencil class is a diagonal scalar stencil, associated with the
// central coefficient of a hex cell. It is defined in the Foam namespace and
// not in the Foam::briscola namespace, just like all other primitives.

namespace Foam
{

class diagStencil;

void writeEntry(Ostream& os, const diagStencil& value);

Istream& operator>>(Istream&, diagStencil&);
Ostream& operator<<(Ostream&, const diagStencil&);

class diagStencil
{

public:

    // The components of this vector space
    scalar v_;

    // Component type
    typedef scalar cmptType;

    // Dimensionality of space
    static const direction dim = 3;

    // Number of components in this vector space
    static const direction nComponents = 1;

    //- Rank of diagStencil is 1
    static const direction rank = 1;

    // Component labeling enumeration (center)
    enum components { C };

    // Static Data Members

    static const char* const typeName;
    static const char* const componentNames[];
    static const diagStencil zero;
    static const diagStencil one;
    static const diagStencil max;
    static const diagStencil min;
    static const diagStencil rootMax;
    static const diagStencil rootMin;

    // Construct null
    inline diagStencil();

    // Construct initialized to zero
    inline diagStencil(const Foam::zero);

    // Construct from a scalar
    inline diagStencil(const scalar&);

    // Construct from Istream
    inline diagStencil(Istream&);

    // Return the number of elements in the diagStencil
    inline static direction size();

    inline const scalar& component(const direction) const;
    inline scalar& component(const direction);

    inline void component(scalar&, const direction) const;
    inline void replace(const direction, const scalar&);

    // Return a diagStencil with all elements being s
    inline static diagStencil uniform(const scalar& s);

    // Access to the center coefficient

    inline scalar& center();
    inline const scalar& center() const;

    // Member Operators

    inline const scalar& operator[](const direction) const;
    inline scalar& operator[](const direction);

    inline void operator+=(const diagStencil&);
    inline void operator-=(const diagStencil&);

    inline void operator=(const scalar);
    inline void operator=(const Foam::zero);

    inline void operator*=(const scalar);
    inline void operator/=(const scalar);

    // IOstream Operators

    friend Istream& operator>>(Istream&, diagStencil&);
    friend Ostream& operator<<(Ostream&, const diagStencil&);
};

// Return a string representation of a diagStencil
word name(const diagStencil&);

template<>
inline bool contiguous<diagStencil>() {return true;}

}

#include "diagStencilI.H"

#endif
