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
    const Cmpt& x0,
    const Cmpt& x1,
    const Cmpt& x2,
    const Cmpt& x3,
    const Cmpt& y0,
    const Cmpt& y1,
    const Cmpt& y2,
    const Cmpt& y3,
    const Cmpt& z0,
    const Cmpt& z1,
    const Cmpt& z2,
    const Cmpt& z3
)
{
    this->v_[X0] = x0;
    this->v_[X1] = x1;
    this->v_[X2] = x2;
    this->v_[X3] = x3;

    this->v_[Y0] = y0;
    this->v_[Y1] = y1;
    this->v_[Y2] = y2;
    this->v_[Y3] = y3;

    this->v_[Z0] = z0;
    this->v_[Z1] = z1;
    this->v_[Z2] = z2;
    this->v_[Z3] = z3;
}

template<class Cmpt>
inline EdgeSpace<Cmpt>::EdgeSpace(Istream& is)
:
    EdgeSpace::csType(is)
{}

}