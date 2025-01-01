#ifndef stencil_H
#define stencil_H

#include "StencilSpace.H"
#include "briscolaUtils.H"
#include "pTransform.H"
#include "contiguous.H"
#include "Row.H"
#include "symmStencil.H"

namespace Foam
{

class stencil
:
    public StencilSpace<stencil, 7>
{

public:

    // Component labeling enumeration
    enum components { C, L, R, B, T, A, F };

    // A stencil is a full stencil type
    typedef stencil fullStencilType;

    // Construct null
    inline stencil();

    // Construct initialised to zero
    inline stencil(const Foam::zero);

    // Construct from cell space
    inline stencil(const CellSpace<stencil,scalar,7>&);

    // Construct from stencil space
    inline stencil(const StencilSpace<stencil,7>&);

    // Construct uniform
    inline stencil(const scalar&);

    // Construct from 7 components
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

    // Construct from Istream
    inline stencil(Istream&);

    // Access to center, left, right, bottom, top, aft and fore element

    inline const scalar& center() const { return this->v_[C]; }
    inline const scalar& left() const { return this->v_[L]; }
    inline const scalar& right() const { return this->v_[R]; }
    inline const scalar& bottom() const { return this->v_[B]; }
    inline const scalar& top() const { return this->v_[T]; }
    inline const scalar& aft() const { return this->v_[A]; }
    inline const scalar& fore() const { return this->v_[F]; }

    inline scalar& center()  { return this->v_[C]; }
    inline scalar& left()  { return this->v_[L]; }
    inline scalar& right()  { return this->v_[R]; }
    inline scalar& bottom()  { return this->v_[B]; }
    inline scalar& top()  { return this->v_[T]; }
    inline scalar& aft()  { return this->v_[A]; }
    inline scalar& fore()  { return this->v_[F]; }

    // Conversion to symmetric stencil
    inline operator symmStencil() const;
};

template<>
inline bool contiguous<stencil>() { return true; }

template<>
inline bool contiguous<Row<stencil,scalar>>() { return true; }

template<>
inline bool contiguous<Row<stencil,vector>>() { return true; }

template<>
inline bool contiguous<Row<stencil,tensor>>() { return true; }

// stencil casts

template<>
class StencilSpaceCast<stencil::ssType,stencil::ssType>
{
public:
    typedef stencil::ssType type;
};

// Stencil partial specialization of pTransform

template<>
class pTransform<stencil>
{
    stencil p_;

public:

    inline pTransform(const labelTensor& T, const stencil& h)
    {
        if (T == briscola::eye)
        {
            p_ = h;
        }
        else
        {
            const labelVector x(T & briscola::unitX);
            const labelVector y(T & briscola::unitY);
            const labelVector z(T & briscola::unitZ);

            p_[0] = h.center();

            p_[briscola::faceNumber(-x)+1] = h.left();
            p_[briscola::faceNumber( x)+1] = h.right();

            p_[briscola::faceNumber(-y)+1] = h.bottom();
            p_[briscola::faceNumber( y)+1] = h.top();

            p_[briscola::faceNumber(-z)+1] = h.aft();
            p_[briscola::faceNumber( z)+1] = h.fore();
        }
    }

    inline operator stencil()
    {
        return p_;
    }
};

// stencil is (quite obviously) a stencil

template<>
inline bool isStencil<stencil>()
{
    return true;
}

}

#include "stencilProducts.H"
#include "stencilI.H"

#endif