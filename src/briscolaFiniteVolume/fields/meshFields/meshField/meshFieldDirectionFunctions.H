// Cell access

inline Type& operator()
(
    const int l,
    const int d,
    const int i,
    const int j,
    const int k
)
{
    return direction(l,d)(i,j,k);
}

inline const Type& operator()
(
    const int l,
    const int d,
    const int i,
    const int j,
    const int k
) const
{
    return direction(l,d)(i,j,k);
}

inline Type& operator()
(
    const int l,
    const int d,
    const labelVector& ijk
)
{
    return direction(l,d)(ijk);
}

inline const Type& operator()
(
    const int l,
    const int d,
    const labelVector& ijk
) const
{
    return direction(l,d)(ijk);
}

// Access to the base block
inline block<Type>& B(const int l, const int d)
{
    return direction(l,d).B();
}

// Const access to the base block
inline const block<Type>& B(const int l, const int d) const
{
    return direction(l,d).B();
}

// Internal cell start and end indices
inline const faceLabel& I(const int l, const int d) const
{
    return direction(l,d).I();
}

// Internal cell block size
inline labelVector N(const int l, const int d) const
{
    return direction(l,d).N();
}
