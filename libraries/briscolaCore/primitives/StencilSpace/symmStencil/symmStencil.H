#ifndef symmStencil_H
#define symmStencil_H

#include "StencilSpace.H"
#include "briscolaUtils.H"
#include "pTransform.H"
#include "contiguous.H"
#include "Row.H"

namespace Foam
{

class stencil;

class symmStencil
:
    public StencilSpace<symmStencil, 4>
{

public:

    // Component labeling enumeration
    enum components { C, L, B, A };

    // A stencil is the full stencil type of symmStencil
    typedef stencil fullStencilType;

    // Construct null
    inline symmStencil();

    // Construct initialised to zero
    inline symmStencil(const Foam::zero);

    // Construct from cell space
    inline symmStencil(const CellSpace<symmStencil,scalar,4>&);

    // Construct from symmStencil space
    inline symmStencil(const StencilSpace<symmStencil,4>&);

    // Construct uniform
    inline symmStencil(const scalar&);

    // Construct from 4 components
    inline symmStencil
    (
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&
    );

    // Construct from Istream
    inline symmStencil(Istream&);

    // Access to center, left, right, bottom, top, aft and fore element

    inline const scalar& center() const { return this->v_[C]; }
    inline const scalar& left() const { return this->v_[L]; }
    inline const scalar& bottom() const { return this->v_[B]; }
    inline const scalar& aft() const { return this->v_[A]; }

    inline scalar& center()  { return this->v_[C]; }
    inline scalar& left()  { return this->v_[L]; }
    inline scalar& bottom()  { return this->v_[B]; }
    inline scalar& aft()  { return this->v_[A]; }
};

template<>
inline bool contiguous<symmStencil>() { return true; }

template<>
inline bool contiguous<Row<symmStencil,scalar>>() { return true; }

template<>
inline bool contiguous<Row<symmStencil,vector>>() { return true; }

template<>
inline bool contiguous<Row<symmStencil,tensor>>() { return true; }

// symmStencil casts

template<>
class StencilSpaceCast<symmStencil::ssType,symmStencil::ssType>
{
public:
    typedef symmStencil::ssType type;
};

// symmStencil is a stencil

template<>
inline bool isStencil<symmStencil>()
{
    return true;
}

}

#include "symmStencilProducts.H"
#include "symmStencilI.H"

#endif
