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
};

}

}

}

#endif
