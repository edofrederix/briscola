// Cell access

inline Type& operator()
(
    const int d,
    const int i,
    const int j,
    const int k
)
{
    return direction(d)(i,j,k);
}

inline const Type& operator()
(
    const int d,
    const int i,
    const int j,
    const int k
) const
{
    return direction(d)(i,j,k);
}

inline Type& operator()
(
    const int d,
    const labelVector& ijk
)
{
    return direction(d)(ijk);
}

inline const Type& operator()
(
    const int d,
    const labelVector& ijk
) const
{
    return direction(d)(ijk);
}

// Access to the base block
inline block<Type>& B(const int d)
{
    return direction(d).B();
}

// Const access to the base block
inline const block<Type>& B(const int d) const
{
    return direction(d).B();
}

// Internal cell start and end indices
inline const faceLabel& I(const int d) const
{
    return direction(d).I();
}

// Internal cell block size
inline labelVector N(const int d) const
{
    return direction(d).N();
}

// Boundary start
inline labelVector boundaryStart
(
    const int d,
    const labelVector& bo
) const
{
    return direction(d).boundaryStart(bo);
}

// Boundary end
inline labelVector boundaryEnd
(
    const int d,
    const labelVector& bo
) const
{
    return direction(d).boundaryEnd(bo);
}

// Boundary size
inline labelVector boundaryN(const int d, const labelVector& bo) const
{
    return direction(d).boundaryN(bo);
}

// Shifted boundary
inline bool shifted(const int d, const labelVector& bo) const
{
    return direction(d).shifted(bo);
}
