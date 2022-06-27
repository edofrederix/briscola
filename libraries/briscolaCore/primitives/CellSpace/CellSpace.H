#ifndef CellSpace_H
#define CellSpace_H

#include "direction.H"
#include "scalar.H"
#include "word.H"
#include "zero.H"

namespace Foam
{

template<class Form, class Cmpt, direction Ncmpts>
class CellSpace;

template<class Form, class Cmpt, direction Ncmpts>
void writeEntry
(
    Ostream& os,
    const CellSpace<Form, Cmpt, Ncmpts>& value
);

template<class Form, class Cmpt, direction Ncmpts>
Istream& operator>>
(
    Istream&,
    CellSpace<Form, Cmpt,Ncmpts>&
);

template<class Form, class Cmpt, direction Ncmpts>
Ostream& operator<<
(
    Ostream&,
    const CellSpace<Form, Cmpt, Ncmpts>&
);

// Class to cast Form1 and Form2 to their equivalent with given Cmpt as
// component type. Cell space casts are necessary for functions that return the
// same form, but with differen component type. Should be specialized by derived
// types.

template<class Form1, class Form2, class Cmpt>
class CellSpaceCast
{
public:
    typedef Form1 type;
};

template<class Form, class Cmpt, direction Ncmpts>
class CellSpace
{
public:

    // The components of this cell space
    Cmpt v_[Ncmpts];

    // CellSpace type
    typedef CellSpace<Form, Cmpt, Ncmpts> csType;

    // Component type
    typedef Cmpt cmptType;

    // Dimensionality of space
    static const direction dim = Cmpt::dim;

    // Number of components in this cell space
    static const direction nComponents = Ncmpts;

    // Static Data Members

    static const char* const typeName;
    static const char* const componentNames[];
    static const Form zero;
    static const Form one;
    static const Form max;
    static const Form min;
    static const Form rootMax;
    static const Form rootMin;

    // Construct null
    inline CellSpace();

    // Construct initialised to zero
    inline CellSpace(const Foam::zero);

    // Construct from Istream
    CellSpace(Istream&);

    // Uniform constructor
    inline CellSpace(const Cmpt&);

    // Construct as copy of another cell space with the same size
    template<class Form2, class Cmpt2>
    inline explicit CellSpace(const CellSpace<Form2, Cmpt2, Ncmpts>&);

    // Return the number of elements in the CellSpace = Ncmpts.
    inline static direction size();

    inline const Cmpt& component(const direction) const;
    inline Cmpt& component(const direction);

    inline void component(Cmpt&, const direction) const;
    inline void replace(const direction, const Cmpt&);

    // Return a CellSpace with all elements = s
    inline static Form uniform(const Cmpt& s);

    // Member Operators

    inline const Cmpt& operator[](const direction) const;
    inline Cmpt& operator[](const direction);

    inline void operator+=(const CellSpace<Form, Cmpt, Ncmpts>&);
    inline void operator-=(const CellSpace<Form, Cmpt, Ncmpts>&);

    inline void operator*=
    (
        const CellSpace
        <
            typename CellSpaceCast<Form,Form,scalar>::type,
            scalar,
            Ncmpts
        >&
    );

    inline void operator/=
    (
        const CellSpace
        <
            typename CellSpaceCast<Form,Form,scalar>::type,
            scalar,
            Ncmpts
        >&
    );

    inline void operator=(const Foam::zero);
    inline void operator=(const Cmpt&);
    inline void operator*=(const scalar);
    inline void operator/=(const scalar);

    inline void operator+=(const Cmpt&);
    inline void operator-=(const Cmpt&);


    // IOstream Operators

    friend Istream& operator>> <Form, Cmpt, Ncmpts>
    (
        Istream&,
        CellSpace<Form, Cmpt, Ncmpts>&
    );

    friend Ostream& operator<< <Form, Cmpt, Ncmpts>
    (
        Ostream&,
        const CellSpace<Form, Cmpt, Ncmpts>&
    );
};

// Return a string representation of a CellSpace
template<class Form, class Cmpt, direction Ncmpts>
word name(const CellSpace<Form, Cmpt, Ncmpts>&);

}

#include "CellSpaceI.H"

#ifdef NoRepository
    #include "CellSpace.C"
#endif

#endif