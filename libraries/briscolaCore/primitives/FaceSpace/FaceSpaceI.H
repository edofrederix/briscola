namespace Foam
{

template<class Cmpt>
inline FaceSpace<Cmpt>::FaceSpace()
{}

template<class Cmpt>
inline FaceSpace<Cmpt>::FaceSpace(const zero)
:
    FaceSpace::csType(Zero)
{}

template<class Cmpt>
inline FaceSpace<Cmpt>::FaceSpace(const Cmpt& v)
:
    FaceSpace::csType(v)
{}

template<class Cmpt>
template<class Cmpt2>
inline FaceSpace<Cmpt>::FaceSpace
(
    const CellSpace<FaceSpace<Cmpt2>, Cmpt2, 6>& vs
)
:
    FaceSpace::csType(vs)
{}

template<class Cmpt>
inline FaceSpace<Cmpt>::FaceSpace
(
    const Cmpt& vl,
    const Cmpt& vr,
    const Cmpt& vb,
    const Cmpt& vt,
    const Cmpt& va,
    const Cmpt& vf
)
{
    this->v_[L] = vl;
    this->v_[R] = vr;
    this->v_[B] = vb;
    this->v_[T] = vt;
    this->v_[A] = va;
    this->v_[F] = vf;
}

template<class Cmpt>
inline FaceSpace<Cmpt>::FaceSpace
(
    const Vector<Cmpt>& lower,
    const Vector<Cmpt>& upper
)
{
    this->v_[L] = lower.x();
    this->v_[R] = upper.x();
    this->v_[B] = lower.y();
    this->v_[T] = upper.y();
    this->v_[A] = lower.z();
    this->v_[F] = upper.z();
}

template<class Cmpt>
inline FaceSpace<Cmpt>::FaceSpace(Istream& is)
:
    FaceSpace::csType(is)
{}

template<>
inline FaceSpace<scalar>::operator stencil() const
{
    stencil r;

    r[stencil::C] = 0;
    r[stencil::L] = this->v_[L];
    r[stencil::R] = this->v_[R];
    r[stencil::B] = this->v_[B];
    r[stencil::T] = this->v_[T];
    r[stencil::A] = this->v_[A];
    r[stencil::F] = this->v_[F];

    return r;
}

template<>
inline FaceSpace<scalar>::operator stencil::ssType() const
{
    stencil::ssType r;

    r[stencil::C] = 0;
    r[stencil::L] = this->v_[L];
    r[stencil::R] = this->v_[R];
    r[stencil::B] = this->v_[B];
    r[stencil::T] = this->v_[T];
    r[stencil::A] = this->v_[A];
    r[stencil::F] = this->v_[F];

    return r;
}

template<>
inline FaceSpace<scalar>::operator stencil::csType() const
{
    stencil::csType r;

    r[stencil::C] = 0;
    r[stencil::L] = this->v_[L];
    r[stencil::R] = this->v_[R];
    r[stencil::B] = this->v_[B];
    r[stencil::T] = this->v_[T];
    r[stencil::A] = this->v_[A];
    r[stencil::F] = this->v_[F];

    return r;
}

// Stencil functions that behave as if FaceSpace is a stencil without central
// coefficient

template<class Cmpt>
inline Cmpt neighborSum
(
    const FaceSpace<Cmpt>& fs
)
{
    return
        fs[FaceSpace<Cmpt>::L]
      + fs[FaceSpace<Cmpt>::R]
      + fs[FaceSpace<Cmpt>::B]
      + fs[FaceSpace<Cmpt>::T]
      + fs[FaceSpace<Cmpt>::A]
      + fs[FaceSpace<Cmpt>::F];
}

template<class Cmpt>
inline Cmpt stencilSum
(
    const FaceSpace<Cmpt>& fs
)
{
    return neighborSum(fs);
}

template<class Cmpt>
inline Cmpt neighborSkewSum
(
    const FaceSpace<Cmpt>& fs
)
{
    return
      - fs[FaceSpace<Cmpt>::L]
      + fs[FaceSpace<Cmpt>::R]
      - fs[FaceSpace<Cmpt>::B]
      + fs[FaceSpace<Cmpt>::T]
      - fs[FaceSpace<Cmpt>::A]
      + fs[FaceSpace<Cmpt>::F];
}

}