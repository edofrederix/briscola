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

template<>
inline scalar neighborSum
(
    const stencil::ssType& ss
)
{
    return
      + ss[stencil::L]
      + ss[stencil::R]
      + ss[stencil::B]
      + ss[stencil::T]
      + ss[stencil::A]
      + ss[stencil::F];
}

template<>
inline scalar neighborSkewSum
(
    const stencil::ssType& ss
)
{
    return
      - ss[stencil::L]
      + ss[stencil::R]
      - ss[stencil::B]
      + ss[stencil::T]
      - ss[stencil::A]
      + ss[stencil::F];
}

}