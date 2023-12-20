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

}