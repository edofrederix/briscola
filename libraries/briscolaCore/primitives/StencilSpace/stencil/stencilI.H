namespace Foam
{

inline stencil::stencil()
{}

inline stencil::stencil(const Foam::zero)
:
    stencil::ssType(Zero)
{}

inline stencil::stencil
(
    const CellSpace<stencil,scalar,7>& cs
)
:
    stencil::ssType(cs)
{}

inline stencil::stencil
(
    const StencilSpace<stencil,7>& ss
)
:
    stencil::ssType(ss)
{}

inline stencil::stencil
(
    const scalar& v
)
:
    stencil::ssType(v)
{}

inline stencil::stencil
(
    const scalar& vc,
    const scalar& vl,
    const scalar& vr,
    const scalar& vb,
    const scalar& vt,
    const scalar& va,
    const scalar& vf
)
{
    this->v_[C] = vc;
    this->v_[L] = vl;
    this->v_[R] = vr;
    this->v_[B] = vb;
    this->v_[T] = vt;
    this->v_[A] = va;
    this->v_[F] = vf;
}

inline stencil::stencil(Istream& is)
:
    stencil::ssType(is)
{}

inline stencil::operator symmStencil() const
{
    return symmStencil
    (
        this->v_[C],
        this->v_[L],
        this->v_[B],
        this->v_[A]
    );
}

// Symmetric stencil operators

// With stencil

inline symmStencil operator+
(
    const stencil::csType& cs1,
    const symmStencil::csType& cs2
)
{
    symmStencil::csType r(cs2);

    r +=
        static_cast<const symmStencil&>
        (
            static_cast<const stencil&>(cs1)
        );

    return r;
}

inline symmStencil operator+
(
    const symmStencil::csType& cs1,
    const stencil::csType& cs2
)
{
    return cs2 + cs1;
}

inline symmStencil operator-
(
    const stencil::csType& cs1,
    const symmStencil::csType& cs2
)
{
    symmStencil::csType r(-cs2);

    r +=
        static_cast<const symmStencil&>
        (
            static_cast<const stencil&>(cs1)
        );

    return r;
}

inline symmStencil operator-
(
    const symmStencil::csType& cs1,
    const stencil::csType& cs2
)
{
    symmStencil::csType r(cs1);

    r -=
        static_cast<const symmStencil&>
        (
            static_cast<const stencil&>(cs2)
        );

    return r;
}

inline symmStencil operator*
(
    const stencil::csType& cs1,
    const symmStencil::csType& cs2
)
{
    symmStencil::csType r(cs2);

    r *=
        static_cast<const symmStencil&>
        (
            static_cast<const stencil&>(cs1)
        );

    return r;
}

inline symmStencil operator*
(
    const symmStencil::csType& cs1,
    const stencil::csType& cs2
)
{
    return cs2 * cs1;
}

inline symmStencil operator/
(
    const stencil::csType& cs1,
    const symmStencil::csType& cs2
)
{
    symmStencil::csType r
    (
        static_cast<const symmStencil&>
        (
            static_cast<const stencil&>(cs1)
        )
    );

    r /= cs2;

    return r;
}

inline symmStencil operator/
(
    const symmStencil::csType& cs1,
    const stencil::csType& cs2
)
{
    symmStencil::csType r(cs1);

    r /=
        static_cast<const symmStencil&>
        (
            static_cast<const stencil&>(cs2)
        );

    return r;
}

}
