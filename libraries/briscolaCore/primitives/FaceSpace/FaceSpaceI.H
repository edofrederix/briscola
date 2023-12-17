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

template<class Cmpt>
inline FaceSpace<Cmpt>::operator LowerFaceSpace<Cmpt>() const
{
    LowerFaceSpace<Cmpt> r;

    r[0] = this->v_[L];
    r[1] = this->v_[B];
    r[2] = this->v_[A];

    return r;
}

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

// Functions with lower face space

// Free functions

template<class Cmpt>
inline LowerFaceSpace<Cmpt> cmptMultiply
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = cmptMultiply(fs[i*2], ls.v_[i]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> cmptMultiply
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = cmptMultiply(ls.v_[i], fs[i*2]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> cmptPow
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = cmptPow(fs[i*2], ls.v_[i]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> cmptPow
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = cmptPow(ls.v_[i], fs[i*2]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> cmptDivide
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = cmptDivide(fs[i*2], ls.v_[i]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> cmptDivide
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = cmptDivide(ls.v_[i], fs[i*2]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> max
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = max(fs[i*2], ls.v_[i]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> max
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = max(ls.v_[i], fs[i*2]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> min
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = min(fs[i*2], ls.v_[i]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> min
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = min(ls.v_[i], fs[i*2]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> minMod
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = minMod(fs[i*2], ls.v_[i]);

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> minMod
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    LowerFaceSpace<Cmpt> r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = minMod(ls.v_[i], fs[i*2]);

    return r;
}

// Operators

template<class Cmpt>
inline LowerFaceSpace<Cmpt> operator+
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    LowerFaceSpace<Cmpt> r(fs);

    r += ls;

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> operator+
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    LowerFaceSpace<Cmpt> r(fs);

    r += ls;

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> operator-
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    LowerFaceSpace<Cmpt> r(fs);

    r -= ls;

    return r;
}

template<class Cmpt>
inline LowerFaceSpace<Cmpt> operator-
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    LowerFaceSpace<Cmpt> r(-fs);

    r += ls;

    return r;
}

template<class Cmpt1, class Cmpt2>
inline typename CellSpaceCast
<
    FaceSpace<Cmpt1>,
    LowerFaceSpace<Cmpt2>,
    typename outerProduct<Cmpt1,Cmpt2>::type
>::type operator*
(
    const FaceSpace<Cmpt1>& fs,
    const LowerFaceSpace<Cmpt2>& ls
)
{
    typename CellSpaceCast
    <
        FaceSpace<Cmpt1>,
        LowerFaceSpace<Cmpt2>,
        typename outerProduct<Cmpt1,Cmpt2>::type
    >::type r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = fs[i*2] * ls.v_[i];

    return r;
}

template<class Cmpt1, class Cmpt2>
inline typename CellSpaceCast
<
    LowerFaceSpace<Cmpt2>,
    FaceSpace<Cmpt1>,
    typename outerProduct<Cmpt2,Cmpt1>::type
>::type operator*
(
    const LowerFaceSpace<Cmpt2>& ls,
    const FaceSpace<Cmpt1>& fs
)
{
    typename CellSpaceCast
    <
        LowerFaceSpace<Cmpt2>,
        FaceSpace<Cmpt1>,
        typename outerProduct<Cmpt2,Cmpt1>::type
    >::type r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = ls.v_[i] * fs[i*2];

    return r;
}

template<class Cmpt1, class Cmpt2>
inline typename CellSpaceCast
<
    FaceSpace<Cmpt1>,
    LowerFaceSpace<Cmpt2>,
    Cmpt1
>::type operator/
(
    const FaceSpace<Cmpt1>& fs,
    const LowerFaceSpace<Cmpt2>& ls
)
{
    typename CellSpaceCast
    <
        FaceSpace<Cmpt1>,
        LowerFaceSpace<Cmpt2>,
        Cmpt1
    >::type r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = fs[i*2]/ls.v_[i];

    return r;
}

template<class Cmpt1, class Cmpt2>
inline typename CellSpaceCast
<
    LowerFaceSpace<Cmpt1>,
    FaceSpace<Cmpt2>,
    Cmpt1
>::type operator/
(
    const LowerFaceSpace<Cmpt1>& ls,
    const FaceSpace<Cmpt2>& fs
)
{
    typename CellSpaceCast
    <
        LowerFaceSpace<Cmpt1>,
        FaceSpace<Cmpt2>,
        Cmpt1
    >::type r;

    for (label i = 0; i < 3; i++)
        r.v_[i] = ls.v_[i]/fs[i*2];

    return r;
}

template<class Cmpt>
inline bool operator==
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    for (label i = 0; i < 3; i++)
        if (fs[i*2] != ls.v_[i])
            return false;

    return true;
}

template<class Cmpt>
inline bool operator==
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    for (label i = 0; i < 3; i++)
        if (fs[i*2] != ls.v_[i])
            return false;

    return true;
}

template<class Cmpt>
inline bool operator!=
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    return !(fs == ls);
}

template<class Cmpt>
inline bool operator!=
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    return !(fs == ls);
}

template<class Cmpt>
inline bool operator>
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    for (label i = 0; i < 3; i++)
        if (fs[i*2] <= ls.v_[i])
            return false;

    return true;
}

template<class Cmpt>
inline bool operator>
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    for (label i = 0; i < 3; i++)
        if (ls.v_[i] <= fs[i*2])
            return false;

    return true;
}

template<class Cmpt>
inline bool operator<
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    for (label i = 0; i < 3; i++)
        if (fs[i*2] >= ls.v_[i])
            return false;

    return true;
}

template<class Cmpt>
inline bool operator<
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    for (label i = 0; i < 3; i++)
        if (ls.v_[i] >= fs[i*2])
            return false;

    return true;
}

template<class Cmpt>
inline bool operator>=
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    return !(fs < ls);
}

template<class Cmpt>
inline bool operator>=
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    return !(ls < fs);
}

template<class Cmpt>
inline bool operator<=
(
    const FaceSpace<Cmpt>& fs,
    const LowerFaceSpace<Cmpt>& ls
)
{
    return !(fs > ls);
}

template<class Cmpt>
inline bool operator<=
(
    const LowerFaceSpace<Cmpt>& ls,
    const FaceSpace<Cmpt>& fs
)
{
    return !(ls > fs);
}

}
