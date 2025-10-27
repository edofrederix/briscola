#ifndef staggered_H
#define staggered_H

#include "meshType.H"
#include "scalar.H"
#include "tensor.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class staggered : public meshType<staggered>
{
public:

    // Vectors are expressed as scalars on staggered meshes
    typedef scalar vectorType;

    // Tensors are expressed as vectors on staggered meshes
    typedef vector tensorType;

    // Dimensionality of staggered space is 3
    static const direction dim = 3;

    // Dimensionality type of staggered is a vector
    typedef vector dimType;

    // Project the given type onto a direction given a base tensor
    template<class Type>
    inline static typename innerProduct<Type,vector>::type project
    (
        const Type& v,
        const direction d,
        const tensor& base
    )
    {
        return (base.vectorComponent(d) & v);
    }
};

}

}

}

#endif
