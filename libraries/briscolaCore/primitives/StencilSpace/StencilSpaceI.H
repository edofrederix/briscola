namespace Foam
{

template<class Form, direction Ncmpts>
inline StencilSpace<Form, Ncmpts>::StencilSpace()
{}

template<class Form, direction Ncmpts>
inline StencilSpace<Form, Ncmpts>::StencilSpace(const zero)
:
    StencilSpace::csType(Zero)
{}

template<class Form, direction Ncmpts>
inline StencilSpace<Form, Ncmpts>::StencilSpace
(
    const scalar& v
)
:
    StencilSpace::csType(v)
{}

template<class Form, direction Ncmpts>
inline StencilSpace<Form, Ncmpts>::StencilSpace
(
    const CellSpace<Form, scalar, Ncmpts>& cs
)
:
    StencilSpace::csType(cs)
{}

template<class Form, direction Ncmpts>
inline StencilSpace<Form, Ncmpts>::StencilSpace(Istream& is)
:
    StencilSpace::csType(is)
{}

}