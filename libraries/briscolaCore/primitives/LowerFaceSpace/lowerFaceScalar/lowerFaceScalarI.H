namespace Foam
{

// Symmetric stencil operators. Arguments need to be defined in terms of cell
// space, to avoid ambiguities. They are then cast to their actual form, where
// the conversion is defined.

inline symmStencil operator+
(
    const symmStencil::csType& cs1,
    const lowerFaceScalar::csType& cs2
)
{
    symmStencil r(cs1);

    r +=
        static_cast<const symmStencil&>
        (
            static_cast<const lowerFaceScalar&>(cs2)
        );

    return r;
}

inline symmStencil operator+
(
    const lowerFaceScalar::csType& cs1,
    const symmStencil::csType& cs2
)
{
    return cs2 + cs1;
}

inline symmStencil operator-
(
    const symmStencil::csType& cs1,
    const lowerFaceScalar::csType& cs2
)
{
    symmStencil::csType r(cs1);

    r -=
        static_cast<const symmStencil&>
        (
            static_cast<const lowerFaceScalar&>(cs2)
        );

    return r;
}

inline symmStencil operator-
(
    const lowerFaceScalar::csType& cs1,
    const symmStencil::csType& cs2
)
{
    symmStencil::csType r(-cs2);

    r +=
        static_cast<const symmStencil&>
        (
            static_cast<const lowerFaceScalar&>(cs1)
        );

    return r;
}

inline symmStencil operator*
(
    const symmStencil::csType& cs1,
    const lowerFaceScalar::csType& cs2
)
{
    symmStencil::csType r(cs1);

    r *=
        static_cast<const symmStencil&>
        (
            static_cast<const lowerFaceScalar&>(cs2)
        );

    return r;
}

inline symmStencil operator*
(
    const lowerFaceScalar::csType& cs2,
    const symmStencil::csType& cs1
)
{
    return cs2 * cs1;
}

inline symmStencil operator/
(
    const lowerFaceScalar::csType& cs1,
    const symmStencil::csType& cs2
)
{
    symmStencil::csType r
    (
        static_cast<const symmStencil&>
        (
            static_cast<const lowerFaceScalar&>(cs1)
        )
    );

    r /= cs2;

    return r;
}

CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pow3)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pow4)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pow5)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pow6)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pow025)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, sqrt)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, cbrt)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, sign)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pos)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pos0)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, neg)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, neg0)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, posPart)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, negPart)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, exp)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, log)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, log10)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, sin)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, cos)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, tan)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, asin)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, acos)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, atan)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, sinh)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, cosh)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, tanh)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, asinh)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, acosh)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, atanh)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, erf)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, erfc)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, lgamma)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, j0)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, j1)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, y0)
CS_UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, y1)

}
