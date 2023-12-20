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

}