namespace Foam
{

inline diagStencil::diagStencil()
{}

inline diagStencil::diagStencil(const Foam::zero)
:
    diagStencil::ssType(Zero)
{}

inline diagStencil::diagStencil
(
    const CellSpace<diagStencil,scalar,1>& cs
)
:
    diagStencil::ssType(cs)
{}

inline diagStencil::diagStencil
(
    const StencilSpace<diagStencil,1>& ss
)
:
    diagStencil::ssType(ss)
{}

inline diagStencil::diagStencil
(
    const scalar& vc
)
{
    this->v_[C] = vc;
}

inline diagStencil::diagStencil(Istream& is)
:
    diagStencil::ssType(is)
{}

inline diagStencil::operator stencil() const
{
    stencil r(Zero);
    r[stencil::C] = this->v_[C];
    return r;
}

template<>
inline scalar neighborSkewSum
(
    const diagStencil::ssType& ss
)
{
    return 0;
}

// Stencil operators. Arguments need to be defined in terms of cell space, to
// avoid ambiguities. They are then cast to their actual form, where the
// conversion is defined.

inline stencil operator+
(
    const stencil::csType& cs1,
    const diagStencil::csType& cs2
)
{
    stencil r(cs1);

    r += static_cast<const diagStencil&>(cs2);

    return r;
}

inline stencil operator+
(
    const diagStencil::csType& cs1,
    const stencil::csType& cs2
)
{
    return cs2 + cs1;
}

inline stencil operator-
(
    const stencil::csType& cs1,
    const diagStencil::csType& cs2
)
{
    stencil::csType r(cs1);

    r -= static_cast<const diagStencil&>(cs2);

    return r;
}

inline stencil operator-
(
    const diagStencil::csType& cs1,
    const stencil::csType& cs2
)
{
    stencil::csType r(-cs2);

    r += static_cast<const diagStencil&>(cs1);

    return r;
}

inline stencil operator*
(
    const stencil::csType& cs1,
    const diagStencil::csType& cs2
)
{
    stencil::csType r(cs1);

    r *= static_cast<const diagStencil&>(cs2);

    return r;
}

inline stencil operator*
(
    const diagStencil::csType& cs2,
    const stencil::csType& cs1
)
{
    return cs2 * cs1;
}

inline stencil operator/
(
    const diagStencil::csType& cs1,
    const stencil::csType& cs2
)
{
    stencil::csType r(static_cast<const diagStencil&>(cs1));

    r /= cs2;

    return r;
}

}