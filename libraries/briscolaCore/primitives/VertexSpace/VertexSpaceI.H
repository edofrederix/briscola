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
    const Cmpt& v0,
    const Cmpt& v1,
    const Cmpt& v2,
    const Cmpt& v3,
    const Cmpt& v4,
    const Cmpt& v5,
    const Cmpt& v6,
    const Cmpt& v7
)
{
    this->v_[V0] = v0;
    this->v_[V1] = v1;
    this->v_[V2] = v2;
    this->v_[V3] = v3;
    this->v_[V4] = v4;
    this->v_[V5] = v5;
    this->v_[V6] = v6;
    this->v_[V7] = v7;
}

template<class Cmpt>
inline VertexSpace<Cmpt>::VertexSpace(Istream& is)
:
    VertexSpace::csType(is)
{}

}