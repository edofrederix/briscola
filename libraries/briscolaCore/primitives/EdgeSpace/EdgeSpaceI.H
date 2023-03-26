namespace Foam
{

template<class Cmpt>
inline EdgeSpace<Cmpt>::EdgeSpace()
{}

template<class Cmpt>
inline EdgeSpace<Cmpt>::EdgeSpace(const zero)
:
    EdgeSpace::csType(Zero)
{}

template<class Cmpt>
inline EdgeSpace<Cmpt>::EdgeSpace(const Cmpt& v)
:
    EdgeSpace::csType(v)
{}

template<class Cmpt>
template<class Cmpt2>
inline EdgeSpace<Cmpt>::EdgeSpace
(
    const CellSpace<EdgeSpace<Cmpt2>, Cmpt2, 12>& vs
)
:
    EdgeSpace::csType(vs)
{}

template<class Cmpt>
inline EdgeSpace<Cmpt>::EdgeSpace
(
    const Cmpt& ba,
    const Cmpt& ta,
    const Cmpt& bf,
    const Cmpt& tf,
    const Cmpt& la,
    const Cmpt& ra,
    const Cmpt& lf,
    const Cmpt& rf,
    const Cmpt& lb,
    const Cmpt& rb,
    const Cmpt& lt,
    const Cmpt& rt
)
{
    this->v_[BA] = ba;
    this->v_[TA] = ta;
    this->v_[BF] = bf;
    this->v_[TF] = tf;
    this->v_[LA] = la;
    this->v_[RA] = ra;
    this->v_[LF] = lf;
    this->v_[RF] = rf;
    this->v_[LB] = lb;
    this->v_[RB] = rb;
    this->v_[LT] = lt;
    this->v_[RT] = rt;
}

template<class Cmpt>
inline EdgeSpace<Cmpt>::EdgeSpace(Istream& is)
:
    EdgeSpace::csType(is)
{}

}