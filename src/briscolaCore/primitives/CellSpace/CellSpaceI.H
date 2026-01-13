#include "error.H"
#include "VectorSpace.H"

namespace Foam
{

template<class Form, class Cmpt, direction Ncmpts>
inline CellSpace<Form, Cmpt, Ncmpts>::CellSpace()
{}

template<class Form, class Cmpt, direction Ncmpts>
inline CellSpace<Form, Cmpt, Ncmpts>::CellSpace(const Foam::zero)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] = Zero;
}

template<class Form, class Cmpt, direction Ncmpts>
inline CellSpace<Form, Cmpt, Ncmpts>::CellSpace
(
    const Cmpt& v
)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] = v;
}

template<class Form, class Cmpt, direction Ncmpts>
template<class Form2, class Cmpt2>
inline CellSpace<Form, Cmpt, Ncmpts>::CellSpace
(
    const CellSpace<Form2, Cmpt2, Ncmpts>& cs
)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] = cs.v_[i];
}

template<class Form, class Cmpt, direction Ncmpts>
inline direction CellSpace<Form, Cmpt, Ncmpts>::size()
{
    return Ncmpts;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form CellSpace<Form, Cmpt, Ncmpts>::uniform(const Cmpt& s)
{
    CellSpace<Form, Cmpt, Ncmpts> v;

    for (label i = 0; i < Ncmpts; i++)
        v.v_[i] = s;

    return v;
}

template<class Form, class Cmpt, direction Ncmpts>
inline const Cmpt& CellSpace<Form, Cmpt, Ncmpts>::operator[]
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= Ncmpts)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

template<class Form, class Cmpt, direction Ncmpts>
inline Cmpt& CellSpace<Form, Cmpt, Ncmpts>::operator[]
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= Ncmpts)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

template<class Form, class Cmpt, direction Ncmpts>
inline void CellSpace<Form, Cmpt, Ncmpts>::operator+=
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] += cs.v_[i];
}

template<class Form, class Cmpt, direction Ncmpts>
inline void CellSpace<Form, Cmpt, Ncmpts>::operator-=
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] -= cs.v_[i];
}

template<class Form, class Cmpt, direction Ncmpts>
inline void CellSpace<Form, Cmpt, Ncmpts>::operator*=
(
    const CellSpace
    <
        typename CellSpaceCast<Form,Form,scalar>::type,
        scalar,
        Ncmpts
    >& cs
)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] *= cs.v_[i];
}

template<class Form, class Cmpt, direction Ncmpts>
inline void CellSpace<Form, Cmpt, Ncmpts>::operator/=
(
    const CellSpace
    <
        typename CellSpaceCast<Form,Form,scalar>::type,
        scalar,
        Ncmpts
    >& cs
)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] /= cs.v_[i];
}

template<class Form, class Cmpt, direction Ncmpts>
inline void CellSpace<Form, Cmpt, Ncmpts>::operator=(const Foam::zero)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] = Zero;
}

template<class Form, class Cmpt, direction Ncmpts>
inline void CellSpace<Form, Cmpt, Ncmpts>::operator=(const Cmpt& v)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] = v;
}

template<class Form, class Cmpt, direction Ncmpts>
inline void CellSpace<Form, Cmpt, Ncmpts>::operator*=
(
    const scalar s
)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] *= s;
}

template<class Form, class Cmpt, direction Ncmpts>
inline void CellSpace<Form, Cmpt, Ncmpts>::operator/=
(
    const scalar s
)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] /= s;
}

template<class Form, class Cmpt, direction Ncmpts>
inline void CellSpace<Form, Cmpt, Ncmpts>::operator+=
(
    const Cmpt& v
)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] += v;
}

template<class Form, class Cmpt, direction Ncmpts>
inline void CellSpace<Form, Cmpt, Ncmpts>::operator-=
(
    const Cmpt& v
)
{
    for (label i = 0; i < Ncmpts; i++)
        v_[i] -= v;
}

template<class Form, class Cmpt, direction Ncmpts>
inline typename CellSpaceCast<Form, Form, scalar>::type magSqr
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    typename CellSpaceCast<Form, Form, scalar>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = magSqr(cs.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline typename CellSpaceCast<Form, Form, scalar>::type mag
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    typename CellSpaceCast<Form, Form, scalar>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = mag(cs.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form normalised
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = normalised(cs.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form cmptMultiply
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cmptMultiply(cs1.v_[i], cs2.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form cmptPow
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cmptPow(cs1.v_[i], cs2.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form cmptDivide
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cmptDivide(cs1.v_[i], cs2.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form stabilise
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs,
    const Cmpt& sml
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = stabilise(cs.v_[i], sml);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type
cmptMax
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cmptMax(cs.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type
cmptMin
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cmptMin(cs.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type
cmptSum
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cmptSum(cs.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type
cmptAv
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cmptAv(cs.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type
cmptProduct
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cmptProduct(cs.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form cmptSqr
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cmptSqr(cs.v_[i]);

    return r;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form cmptMag
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cmptMag(cs.v_[i]);

    return r;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form max
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = max(cs1.v_[i], cs2.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form max
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const Cmpt& s2
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = max(cs1.v_[i], s2);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form max
(
    const Cmpt& s1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = max(s1, cs2.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form min
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = min(cs1.v_[i], cs2.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form min
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const Cmpt& s2
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = min(cs1.v_[i], s2);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form min
(
    const Cmpt& s1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = min(s1, cs2.v_[i]);

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form minMod
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = minMod(cs1.v_[i], cs2.v_[i]);

    return r;
}

template
<
    class Form1, class Cmpt1,
    class Form2, class Cmpt2,
    direction Ncmpts
>
inline typename innerProduct
<
    CellSpace<Form1, Cmpt1, Ncmpts>,
    CellSpace<Form2, Cmpt2, Ncmpts>
>::type dot
(
    const CellSpace<Form1, Cmpt1, Ncmpts>& t1,
    const CellSpace<Form2, Cmpt2, Ncmpts>& t2
)
{
    typename innerProduct
    <
        CellSpace<Form1, Cmpt1, Ncmpts>,
        CellSpace<Form2, Cmpt2, Ncmpts>
    >::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = dot(t1.v_[i], t2.v_[i]);

    return r;
}

// Operators

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator-
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = -cs.v_[i];

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator+
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    CellSpace<Form, Cmpt, Ncmpts> r(cs1);

    r += cs2;

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator+
(
    const Cmpt& s,
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    CellSpace<Form, Cmpt, Ncmpts> r(cs);

    r += s;

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator+
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs,
    const Cmpt& s
)
{
    CellSpace<Form, Cmpt, Ncmpts> r(cs);

    r += s;

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator-
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    CellSpace<Form, Cmpt, Ncmpts> r(cs1);

    r -= cs2;

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator-
(
    const Cmpt& s,
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    CellSpace<Form, Cmpt, Ncmpts> r(-cs);

    r += s;

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator-
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs,
    const Cmpt& s
)
{
    CellSpace<Form, Cmpt, Ncmpts> r(cs);

    r -= s;

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator*
(
    scalar s,
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    CellSpace<Form, Cmpt, Ncmpts> r(cs);

    r *= s;

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator*
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs,
    scalar s
)
{
    return s*cs;
}

template<class Form, direction Ncmpts>
inline Form sqr
(
    const CellSpace<Form, scalar, Ncmpts>& cs
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = sqr(cs.v_[i]);

    return r;
}

template<class Form, direction Ncmpts>
inline Form pow
(
    const CellSpace<Form, scalar, Ncmpts>& cs,
    scalar p
)
{
    Form r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = pow(cs.v_[i], p);

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts
>
inline typename CellSpaceCast
<
    Form1,
    Form2,
    typename outerProduct<Cmpt1,Cmpt2>::type
>::type operator*
(
    const CellSpace<Form1, Cmpt1, Ncmpts>& cs1,
    const CellSpace<Form2, Cmpt2, Ncmpts>& cs2
)
{
    typename CellSpaceCast
    <
        Form1,
        Form2,
        typename outerProduct<Cmpt1,Cmpt2>::type
    >::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cs1.v_[i] * cs2.v_[i];

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts1,
    direction Ncmpts2
>
inline typename CellSpaceCast
<
    Form1,
    Form1,
    typename outerProduct<Cmpt1,Form2>::type
>::type operator*
(
    const CellSpace<Form1, Cmpt1, Ncmpts1>& cs1,
    const VectorSpace<Form2, Cmpt2, Ncmpts2>& v2
)
{
    typename CellSpaceCast
    <
        Form1,
        Form1,
        typename outerProduct<Cmpt1,Form2>::type
    >::type r;

    for (label i = 0; i < Ncmpts1; i++)
        r.v_[i] = cs1.v_[i] * static_cast<const Form2&>(v2);

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts1,
    direction Ncmpts2
>
inline typename CellSpaceCast
<
    Form2,
    Form2,
    typename outerProduct<Form1,Cmpt2>::type
>::type operator*
(
    const VectorSpace<Form1, Cmpt1, Ncmpts1>& v1,
    const CellSpace<Form2, Cmpt2, Ncmpts2>& cs2
)
{
    typename CellSpaceCast
    <
        Form2,
        Form2,
        typename outerProduct<Form1,Cmpt2>::type
    >::type r;

    for (label i = 0; i < Ncmpts2; i++)
        r.v_[i] = static_cast<const Form1&>(v1) * cs2.v_[i];

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator/
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs,
    const scalar s
)
{
    CellSpace<Form, Cmpt, Ncmpts> r(cs);

    r /= s;

    return r;
}

template<class Form, direction Ncmpts>
inline Form operator/
(
    const scalar s,
    const CellSpace<Form, scalar, Ncmpts>& cs
)
{
    CellSpace<Form, scalar, Ncmpts> r(s);

    r /= cs;

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts1,
    direction Ncmpts2
>
inline typename CellSpaceCast<Form2, Form2, Form1>::type operator/
(
    const VectorSpace<Form1, Cmpt1, Ncmpts1>& v,
    const CellSpace<Form2, Cmpt2, Ncmpts2>& cs
)
{
    typename CellSpaceCast<Form2, Form2, Form1>::type r;

    for (label i = 0; i < Ncmpts2; i++)
        r.v_[i] = v/cs.v_[i];

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts
>
inline typename CellSpaceCast<Form1, Form2, Cmpt1>::type operator/
(
    const CellSpace<Form1, Cmpt1, Ncmpts>& cs1,
    const CellSpace<Form2, Cmpt2, Ncmpts>& cs2
)
{
    typename CellSpaceCast<Form1, Form2, Cmpt1>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cs1.v_[i]/cs2.v_[i];

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts
>
inline typename CellSpaceCast
<
    Form1,
    Form2,
    typename innerProduct<Cmpt1,Cmpt2>::type
>::type operator&
(
    const CellSpace<Form1, Cmpt1, Ncmpts>& cs1,
    const CellSpace<Form2, Cmpt2, Ncmpts>& cs2
)
{
    typename CellSpaceCast
    <
        Form1,
        Form2,
        typename innerProduct<Cmpt1,Cmpt2>::type
    >::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cs1.v_[i] & cs2.v_[i];

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts1,
    direction Ncmpts2
>
inline typename CellSpaceCast
<
    Form1,
    Form1,
    typename innerProduct<Cmpt1,Form2>::type
>::type operator&
(
    const CellSpace<Form1, Cmpt1, Ncmpts1>& cs1,
    const VectorSpace<Form2, Cmpt2, Ncmpts2>& v2
)
{
    typename CellSpaceCast
    <
        Form1,
        Form1,
        typename innerProduct<Cmpt1,Form2>::type
    >::type r;

    for (label i = 0; i < Ncmpts1; i++)
        r.v_[i] = cs1.v_[i] & static_cast<const Form2&>(v2);

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts1,
    direction Ncmpts2
>
inline typename CellSpaceCast
<
    Form2,
    Form2,
    typename innerProduct<Form1,Cmpt2>::type
>::type
operator&
(
    const VectorSpace<Form1, Cmpt1, Ncmpts1>& v1,
    const CellSpace<Form2, Cmpt2, Ncmpts2>& cs2
)
{
    typename CellSpaceCast
    <
        Form2,
        Form2,
        typename innerProduct<Form1,Cmpt2>::type
    >::type r;

    for (label i = 0; i < Ncmpts2; i++)
        r.v_[i] = static_cast<const Form1&>(v1) & cs2.v_[i];

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline typename CellSpaceCast
<
    Form,
    Form,
    typename Form::cmptCmpt
>::type
operator&&
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cs1.v_[i] && cs2.v_[i];

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type
operator&&
(
    const Cmpt& v1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = v1 && cs2.v_[i];

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type
operator&&
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const Cmpt& v2
)
{
    typename CellSpaceCast<Form, Form, typename Form::cmptCmpt>::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cs1.v_[i] && v2;

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts
>
inline typename CellSpaceCast
<
    Form1,
    Form2,
    typename crossProduct<Cmpt1,Cmpt2>::type
>::type operator^
(
    const CellSpace<Form1, Cmpt1, Ncmpts>& cs1,
    const CellSpace<Form2, Cmpt2, Ncmpts>& cs2
)
{
    typename CellSpaceCast
    <
        Form1,
        Form2,
        typename crossProduct<Cmpt1,Cmpt2>::type
    >::type r;

    for (label i = 0; i < Ncmpts; i++)
        r.v_[i] = cs1.v_[i] ^ cs2.v_[i];

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts1,
    direction Ncmpts2
>
inline typename CellSpaceCast
<
    Form1,
    Form1,
    typename crossProduct<Cmpt1,Form2>::type
>::type
operator^
(
    const CellSpace<Form1, Cmpt1, Ncmpts1>& cs1,
    const VectorSpace<Form2, Cmpt2, Ncmpts2>& v2
)
{
    typename CellSpaceCast
    <
        Form1,
        Form1,
        typename crossProduct<Cmpt1,Form2>::type
    >::type r;

    for (label i = 0; i < Ncmpts1; i++)
        r.v_[i] = cs1.v_[i] ^ static_cast<const Form2&>(v2);

    return r;
}

template
<
    class Form1,
    class Form2,
    class Cmpt1,
    class Cmpt2,
    direction Ncmpts1,
    direction Ncmpts2
>
inline typename CellSpaceCast
<
    Form2,
    Form2,
    typename crossProduct<Form1,Cmpt2>::type
>::type
operator^
(
    const VectorSpace<Form1, Cmpt1, Ncmpts1>& v1,
    const CellSpace<Form2, Cmpt2, Ncmpts2>& cs2
)
{
    typename CellSpaceCast
    <
        Form2,
        Form2,
        typename crossProduct<Form1,Cmpt2>::type
    >::type r;

    for (label i = 0; i < Ncmpts2; i++)
        r.v_[i] = static_cast<const Form1&>(v1) ^ cs2.v_[i];

    return r;
}

template<class Form, class Cmpt, direction Ncmpts>
inline bool operator==
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    for (label i = 0; i < Ncmpts; i++)
        if (cs1.v_[i] != cs2.v_[i])
            return false;

    return true;
}

template<class Form, class Cmpt, direction Ncmpts>
inline bool operator!=
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    return !(cs1 == cs2);
}

template<class Form, class Cmpt, direction Ncmpts>
inline bool operator>
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    for (label i = 0; i < Ncmpts; i++)
        if (cs1.v_[i] <= cs2.v_[i])
            return false;

    return true;
}


template<class Form, class Cmpt, direction Ncmpts>
inline bool operator<
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    for (label i = 0; i < Ncmpts; i++)
        if (cs1.v_[i] >= cs2.v_[i])
            return false;

    return true;
}

template<class Form, class Cmpt, direction Ncmpts>
inline bool operator>=
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    return !(cs1 < cs2);
}

template<class Form, class Cmpt, direction Ncmpts>
inline bool operator<=
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs1,
    const CellSpace<Form, Cmpt, Ncmpts>& cs2
)
{
    return !(cs1 > cs2);
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form trimPrecision(CellSpace<Form, Cmpt, Ncmpts> cs)
{
    for (direction i = 0; i < Ncmpts; i++)
        cs.v_[i] = trimPrecision(cs.v_[i]);

    return cs;
}

}