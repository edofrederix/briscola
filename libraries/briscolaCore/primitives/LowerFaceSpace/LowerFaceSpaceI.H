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
inline LowerFaceSpace<Cmpt>::LowerFaceSpace(const Cmpt& v)
:
    LowerFaceSpace::csType(v)
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

}
