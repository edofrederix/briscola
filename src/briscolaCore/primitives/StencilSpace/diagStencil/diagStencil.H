#ifndef diagStencil_H
#define diagStencil_H

#include "StencilSpace.H"
#include "utils.H"
#include "pTransform.H"
#include "contiguous.H"
#include "Row.H"
#include "stencil.H"

namespace Foam
{

class diagStencil
:
    public StencilSpace<diagStencil, 1>
{

public:

    // Component labeling enumeration
    enum components { C };

    // Construct null
    inline diagStencil();

    // Construct initialised to zero
    inline diagStencil(const Foam::zero);

    // Construct from cell space
    inline diagStencil(const CellSpace<diagStencil,scalar,1>&);

    // Construct from stencil space
    inline diagStencil(const StencilSpace<diagStencil,1>&);

    // Construct from 1 component
    inline diagStencil
    (
        const scalar&
    );

    // Construct from Istream
    inline diagStencil(Istream&);

    // Access to center
    inline const scalar& center() const { return this->v_[C]; }
    inline scalar& center()  { return this->v_[C]; }

    // Conversion to stencil
    inline operator stencil() const;
};

template<>
inline bool contiguous<diagStencil>() { return true; }

template<>
inline bool contiguous<Row<diagStencil,scalar>>() { return true; }

template<>
inline bool contiguous<Row<diagStencil,vector>>() { return true; }

template<>
inline bool contiguous<Row<diagStencil,tensor>>() { return true; }

// diagStencil casts

template<>
class StencilSpaceCast<diagStencil::ssType,diagStencil::ssType>
{
public:
    typedef diagStencil::ssType type;
};

template<>
class StencilSpaceCast<stencil::ssType,diagStencil::ssType>
{
public:
    typedef stencil::ssType type;
};

template<>
class StencilSpaceCast<diagStencil::ssType,stencil::ssType>
{
public:
    typedef stencil::ssType type;
};

// diagStencil partial specialization of pTransform

template<>
class pTransform<diagStencil>
{
    diagStencil p_;

public:

    inline pTransform(const labelTensor&, const diagStencil& h)
    {
        p_ = h;
    }

    inline operator diagStencil()
    {
        return p_;
    }
};

// diagStencil is a stencil

template<>
inline bool isStencil<diagStencil>()
{
    return true;
}

}

#include "diagStencilProducts.H"
#include "diagStencilI.H"

#endif
