namespace Foam
{

template<class Cmpt>
inline LowerFaceSpace<Cmpt>::LowerFaceSpace()
{}

template<class Cmpt>
inline LowerFaceSpace<Cmpt>::LowerFaceSpace(const zero)
:
    LowerFaceSpace::csType(Zero)
{}

template<class Cmpt>
template<class Cmpt2>
inline LowerFaceSpace<Cmpt>::LowerFaceSpace
(
    const CellSpace<LowerFaceSpace<Cmpt2>, Cmpt2, 3>& vs
)
:
    LowerFaceSpace::csType(vs)
{}

template<class Cmpt>
inline LowerFaceSpace<Cmpt>::LowerFaceSpace
(
    const Cmpt& vl,
    const Cmpt& vb,
    const Cmpt& va
)
{
    this->v_[L] = vl;
    this->v_[B] = vb;
    this->v_[A] = va;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt>::LowerFaceSpace
(
    const Vector<Cmpt>& lower
)
{
    this->v_[L] = lower.x();
    this->v_[B] = lower.y();
    this->v_[A] = lower.z();
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt>::LowerFaceSpace(Istream& is)
:
    LowerFaceSpace::csType(is)
{}

template<>
inline LowerFaceSpace<scalar>::operator symmStencil() const
{
    symmStencil r;

    r[symmStencil::C] = 0;
    r[symmStencil::L] = this->v_[L];
    r[symmStencil::B] = this->v_[B];
    r[symmStencil::A] = this->v_[A];

    return r;
}

// Symmetric stencil functions that behave as if LowerFaceSpace is a symmetric
// stencil without central coefficient

template<class Cmpt>
inline Cmpt neighborSum
(
    const LowerFaceSpace<Cmpt>& fs
)
{
    return
        2.0*fs[LowerFaceSpace<Cmpt>::L]
      + 2.0*fs[LowerFaceSpace<Cmpt>::B]
      + 2.0*fs[LowerFaceSpace<Cmpt>::A];
}

template<class Cmpt>
inline Cmpt stencilSum
(
    const LowerFaceSpace<Cmpt>& fs
)
{
    return neighborSum(fs);
}

template<class Cmpt>
inline Cmpt neighborSkewSum
(
    const LowerFaceSpace<Cmpt>& fs
)
{
    return Zero;
}

}
