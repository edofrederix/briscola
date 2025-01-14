// Cell access

inline Type& operator()
(
    const int i,
    const int j,
    const int k
)
{
    return direction()(i,j,k);
}

inline const Type& operator()
(
    const int i,
    const int j,
    const int k
) const
{
    return direction()(i,j,k);
}

inline Type& operator()
(
    const labelVector& ijk
)
{
    return direction()(ijk);
}

inline const Type& operator()
(
    const labelVector& ijk
) const
{
    return direction()(ijk);
}

// Access to the base block
inline block<Type>& B()
{
    return direction().B();
}

// Const access to the base block
inline const block<Type>& B() const
{
    return direction().B();
}

// Internal cell start and end indices
inline const faceLabel& I() const
{
    return direction().I();
}

// Internal cell block size
inline labelVector N() const
{
    return direction().N();
}
