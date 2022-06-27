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

template<class Form, direction Ncmpts>
inline scalar neighborSum
(
    const StencilSpace<Form, Ncmpts>& ss
)
{
    scalar r(0.0);

    for (label i = 1; i < Ncmpts; i++)
        r += ss.v_[i];

    return r;
}

template<class Form, direction Ncmpts>
inline scalar stencilSum
(
    const StencilSpace<Form, Ncmpts>& ss
)
{
    scalar r(ss.v_[0]);

    for (label i = 1; i < Ncmpts; i++)
        r += ss.v_[i];

    return r;
}

template<class Form, direction Ncmpts>
inline scalar neighborSkewSum
(
    const StencilSpace<Form, Ncmpts>& ss
)
{
    // Should be defined by each stencil

    NotImplemented;

    return 0.0;
}

}