namespace Foam
{

// Helper function to transfer row data, allowing for partial specialization

template<class SType, class Type>
inline void transferRowData
(
    Row<SType,Type>& row,
    const SType& s,
    const Type& b
)
{
    for (int i = 0; i < SType::nComponents; i++)
        row[i] = s[i];

    for (int i = 0; i < pTraits<Type>::nComponents; i++)
        row[i+SType::nComponents] = b[i];
}

// Scalar specialization

template<class SType>
inline void transferRowData
(
    Row<SType,scalar>& row,
    const SType& s,
    const scalar& b
)
{
    for (int i = 0; i < SType::nComponents; i++)
        row[i] = s[i];

    row[SType::nComponents] = b;
}

// Helper function to return the source, allowing for partial specialization

template<class SType, class Type>
inline Type sourceFromRow(const Row<SType,Type>& row)
{
    Type b;
    for (int i = 0; i < pTraits<Type>::nComponents; i++)
        b[i] = row[i+SType::nComponents];

    return b;
}

// Scalar specialization

template<class SType>
inline scalar sourceFromRow(const Row<SType,scalar>& row)
{
    return row[SType::nComponents];
}

// Row class

template<class SType, class Type>
inline Row<SType,Type>::Row()
:
    listType()
{}

template<class SType, class Type>
inline Row<SType,Type>::Row(const SType& s, const Type& b)
:
    listType()
{
    transferRowData(*this,s,b);
}

template<class SType, class Type>
inline SType Row<SType,Type>::stencil() const
{
    SType s;
    for (int i = 0; i < SType::nComponents; i++)
        s[i] = this->operator[](i);

    return s;
}

template<class SType, class Type>
inline Type Row<SType,Type>::source() const
{
    return sourceFromRow(*this);
}

}
