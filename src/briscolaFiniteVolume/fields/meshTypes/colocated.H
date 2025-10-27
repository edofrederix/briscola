#ifndef colocated_H
#define colocated_H

#include "meshType.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class colocated : public meshType<colocated>
{
public:

    // Vectors are (trivially) expressed as vectors on colocated meshes
    typedef vector vectorType;

    // Tensors are (trivially) expressed as tensors on colocated meshes
    typedef tensor tensorType;

    // Dimensionality of colocated space is 1
    static const direction dim = 1;

    // Dimensionality type of colocated is a scalar
    typedef scalar dimType;

    // Project the given type onto a direction given a base tensor (this is a
    // placeholder because for colocated meshes there's nothing to project)
    template<class Type>
    inline static Type project
    (
        const Type& v,
        const direction d,
        const tensor& base
    )
    {
        return v;
    }
};

}

}

}

#endif
