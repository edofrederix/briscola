#ifndef pTransform_H
#define pTransform_H

#include "labelTensor.H"

namespace Foam
{

// Default primitive transformation class, which returns an unmodified primitive

template<class Type>
class pTransform
:
    public Type
{
public:

    // Construct from transformation tensor and primitive
    inline pTransform(const labelTensor&, const Type& p)
    :
        Type(p)
    {}

    // Return transformed primitive
    inline operator const Type& () const
    {
        return *this;
    }
};

// Label specialization

template<>
class pTransform<label>
{
    label p_;

public:

    // Construct from transformation tensor and primitive
    inline pTransform(const labelTensor&, const label& p)
    :
        p_(p)
    {}

    // Return transformed primitive
    inline operator label()
    {
        return p_;
    }
};

// Scalar specialization

template<>
class pTransform<scalar>
{
    scalar p_;

public:

    // Construct from transformation tensor and primitive
    inline pTransform(const labelTensor&, const scalar& p)
    :
        p_(p)
    {}

    // Return transformed primitive
    inline operator scalar()
    {
        return p_;
    }
};

}

#endif
