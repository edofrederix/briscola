namespace Foam
{

// Stencil operators. Arguments need to be defined in terms of cell space, to
// avoid ambiguities. They are then cast to their actual form, where the
// conversion is defined.

inline stencil operator+
(
    const stencil::csType& cs1,
    const faceScalar::csType& cs2
)
{
    stencil r(cs1);

    r +=
        static_cast<const stencil&>
        (
            static_cast<const faceScalar&>(cs2)
        );

    return r;
}

inline stencil operator+
(
    const faceScalar::csType& cs1,
    const stencil::csType& cs2
)
{
    return cs2 + cs1;
}

inline stencil operator-
(
    const stencil::csType& cs1,
    const faceScalar::csType& cs2
)
{
    stencil::csType r(cs1);

    r -=
        static_cast<const stencil&>
        (
            static_cast<const faceScalar&>(cs2)
        );

    return r;
}

inline stencil operator-
(
    const faceScalar::csType& cs1,
    const stencil::csType& cs2
)
{
    stencil::csType r(-cs2);

    r +=
        static_cast<const stencil&>
        (
            static_cast<const faceScalar&>(cs1)
        );

    return r;
}

inline stencil operator*
(
    const stencil::csType& cs1,
    const faceScalar::csType& cs2
)
{
    stencil::csType r(cs1);

    r *=
        static_cast<const stencil&>
        (
            static_cast<const faceScalar&>(cs2)
        );

    return r;
}

inline stencil operator*
(
    const faceScalar::csType& cs2,
    const stencil::csType& cs1
)
{
    return cs2 * cs1;
}

inline stencil operator/
(
    const faceScalar::csType& cs1,
    const stencil::csType& cs2
)
{
    stencil::csType r
    (
        static_cast<const stencil&>
        (
            static_cast<const faceScalar&>(cs1)
        )
    );

    r /= cs2;

    return r;
}

CS_UNARY_FUNCTION(faceScalar, faceScalar, pow3)
CS_UNARY_FUNCTION(faceScalar, faceScalar, pow4)
CS_UNARY_FUNCTION(faceScalar, faceScalar, pow5)
CS_UNARY_FUNCTION(faceScalar, faceScalar, pow6)
CS_UNARY_FUNCTION(faceScalar, faceScalar, pow025)
CS_UNARY_FUNCTION(faceScalar, faceScalar, sqrt)
CS_UNARY_FUNCTION(faceScalar, faceScalar, cbrt)
CS_UNARY_FUNCTION(faceScalar, faceScalar, sign)
CS_UNARY_FUNCTION(faceScalar, faceScalar, pos)
CS_UNARY_FUNCTION(faceScalar, faceScalar, pos0)
CS_UNARY_FUNCTION(faceScalar, faceScalar, neg)
CS_UNARY_FUNCTION(faceScalar, faceScalar, neg0)
CS_UNARY_FUNCTION(faceScalar, faceScalar, posPart)
CS_UNARY_FUNCTION(faceScalar, faceScalar, negPart)
CS_UNARY_FUNCTION(faceScalar, faceScalar, exp)
CS_UNARY_FUNCTION(faceScalar, faceScalar, log)
CS_UNARY_FUNCTION(faceScalar, faceScalar, log10)
CS_UNARY_FUNCTION(faceScalar, faceScalar, sin)
CS_UNARY_FUNCTION(faceScalar, faceScalar, cos)
CS_UNARY_FUNCTION(faceScalar, faceScalar, tan)
CS_UNARY_FUNCTION(faceScalar, faceScalar, asin)
CS_UNARY_FUNCTION(faceScalar, faceScalar, acos)
CS_UNARY_FUNCTION(faceScalar, faceScalar, atan)
CS_UNARY_FUNCTION(faceScalar, faceScalar, sinh)
CS_UNARY_FUNCTION(faceScalar, faceScalar, cosh)
CS_UNARY_FUNCTION(faceScalar, faceScalar, tanh)
CS_UNARY_FUNCTION(faceScalar, faceScalar, asinh)
CS_UNARY_FUNCTION(faceScalar, faceScalar, acosh)
CS_UNARY_FUNCTION(faceScalar, faceScalar, atanh)
CS_UNARY_FUNCTION(faceScalar, faceScalar, erf)
CS_UNARY_FUNCTION(faceScalar, faceScalar, erfc)
CS_UNARY_FUNCTION(faceScalar, faceScalar, lgamma)
CS_UNARY_FUNCTION(faceScalar, faceScalar, j0)
CS_UNARY_FUNCTION(faceScalar, faceScalar, j1)
CS_UNARY_FUNCTION(faceScalar, faceScalar, y0)
CS_UNARY_FUNCTION(faceScalar, faceScalar, y1)

}
