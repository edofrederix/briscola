namespace Foam
{

inline symmStencil::symmStencil()
{}

inline symmStencil::symmStencil(const Foam::zero)
:
    symmStencil::ssType(Zero)
{}

inline symmStencil::symmStencil
(
    const CellSpace<symmStencil,scalar,4>& cs
)
:
    symmStencil::ssType(cs)
{}

inline symmStencil::symmStencil
(
    const StencilSpace<symmStencil,4>& ss
)
:
    symmStencil::ssType(ss)
{}

inline symmStencil::symmStencil
(
    const scalar& v
)
:
    symmStencil::ssType(v)
{}

inline symmStencil::symmStencil
(
    const scalar& vc,
    const scalar& vl,
    const scalar& vb,
    const scalar& va
)
{
    this->v_[C] = vc;
    this->v_[L] = vl;
    this->v_[B] = vb;
    this->v_[A] = va;
}

inline symmStencil::symmStencil(Istream& is)
:
    symmStencil::ssType(is)
{}

inline symmStencil::operator stencil() const
{
    stencil r(Zero);
    r[stencil::C] = this->v_[C];
    r[stencil::L] = this->v_[L];
    r[stencil::B] = this->v_[B];
    r[stencil::A] = this->v_[A];
    return r;
}

template<>
inline scalar neighborSum
(
    const symmStencil::ssType& ss
)
{
    return
      + 2.0*ss[symmStencil::L]
      + 2.0*ss[symmStencil::B]
      + 2.0*ss[symmStencil::A];
}

template<>
inline scalar stencilSum
(
    const symmStencil::ssType& ss
)
{
    return
            ss[symmStencil::C]
      + 2.0*ss[symmStencil::L]
      + 2.0*ss[symmStencil::B]
      + 2.0*ss[symmStencil::A];
}

template<>
inline scalar neighborSkewSum
(
    const symmStencil::ssType& ss
)
{
    return 0.0;
}

// Stencil operators

// With stencil

inline stencil operator+
(
    const stencil::csType& cs1,
    const symmStencil::csType& cs2
)
{
    stencil r(cs1);

    r += static_cast<const symmStencil&>(cs2);

    return r;
}

inline stencil operator+
(
    const symmStencil::csType& cs1,
    const stencil::csType& cs2
)
{
    return cs2 + cs1;
}

inline stencil operator-
(
    const stencil::csType& cs1,
    const symmStencil::csType& cs2
)
{
    stencil::csType r(cs1);

    r -= static_cast<const symmStencil&>(cs2);

    return r;
}

inline stencil operator-
(
    const symmStencil::csType& cs1,
    const stencil::csType& cs2
)
{
    stencil::csType r(-cs2);

    r += static_cast<const symmStencil&>(cs1);

    return r;
}

inline stencil operator*
(
    const stencil::csType& cs1,
    const symmStencil::csType& cs2
)
{
    stencil::csType r(cs1);

    r *= static_cast<const symmStencil&>(cs2);

    return r;
}

inline stencil operator*
(
    const symmStencil::csType& cs2,
    const stencil::csType& cs1
)
{
    return cs2 * cs1;
}

inline stencil operator/
(
    const symmStencil::csType& cs1,
    const stencil::csType& cs2
)
{
    stencil::csType r(static_cast<const symmStencil&>(cs1));

    r /= cs2;

    return r;
}

}