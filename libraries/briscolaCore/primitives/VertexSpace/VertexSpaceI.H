namespace Foam
{

template<class Cmpt>
inline VertexSpace<Cmpt>::VertexSpace()
{}

template<class Cmpt>
inline VertexSpace<Cmpt>::VertexSpace(const zero)
:
    VertexSpace::csType(Zero)
{}

template<class Cmpt>
inline VertexSpace<Cmpt>::VertexSpace(const Cmpt& v)
:
    VertexSpace::csType(v)
{}

template<class Cmpt>
template<class Cmpt2>
inline VertexSpace<Cmpt>::VertexSpace
(
    const CellSpace<VertexSpace<Cmpt2>, Cmpt2, 8>& vs
)
:
    VertexSpace::csType(vs)
{}

template<class Cmpt>
inline VertexSpace<Cmpt>::VertexSpace
(
    const Cmpt& lba,
    const Cmpt& rba,
    const Cmpt& lta,
    const Cmpt& rta,
    const Cmpt& lbf,
    const Cmpt& rbf,
    const Cmpt& ltf,
    const Cmpt& rtf
)
{
    this->v_[LBA] = lba;
    this->v_[RBA] = rba;
    this->v_[LTA] = lta;
    this->v_[RTA] = rta;
    this->v_[LBF] = lbf;
    this->v_[RBF] = rbf;
    this->v_[LTF] = ltf;
    this->v_[RTF] = rtf;
}

template<class Cmpt>
inline VertexSpace<Cmpt>::VertexSpace(Istream& is)
:
    VertexSpace::csType(is)
{}

}