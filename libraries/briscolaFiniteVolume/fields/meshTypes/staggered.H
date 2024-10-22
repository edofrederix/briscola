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
};

}

}

}

#endif
