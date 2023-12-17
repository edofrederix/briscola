namespace Foam
{

inline symmStencil::symmStencil()
{}

inline symmStencil::symmStencil(const Foam::zero)
:
    symmStencil::ssType(Zero)
{}

inline symmStencil::symmStencil
(
    const CellSpace<symmStencil,scalar,4>& cs
)
:
    symmStencil::ssType(cs)
{}

inline symmStencil::symmStencil
(
    const StencilSpace<symmStencil,4>& ss
)
:
    symmStencil::ssType(ss)
{}

inline symmStencil::symmStencil
(
    const scalar& v
)
:
    symmStencil::ssType(v)
{}

inline symmStencil::symmStencil
(
    const scalar& vc,
    const scalar& vl,
    const scalar& vb,
    const scalar& va
)
{
    this->v_[C] = vc;
    this->v_[L] = vl;
    this->v_[B] = vb;
    this->v_[A] = va;
}

inline symmStencil::symmStencil(Istream& is)
:
    symmStencil::ssType(is)
{}

}