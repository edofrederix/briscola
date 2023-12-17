#ifndef meshDirectionStencilOperations_H
#define meshDirectionStencilOperations_H

#include "stencil.H"
#include "meshDirectionStencilProductFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Full stencil conversion

template<class MeshType>
inline stencil fullStencil
(
    const meshDirection<stencil,MeshType>& st,
    const label& i,
    const label& j,
    const label& k
)
{
    return st(i,j,k);
}

template<class MeshType>
inline stencil fullStencil
(
    const meshDirection<stencil,MeshType>& st,
    const labelVector ijk
)
{
    return st(ijk);
}

// Stencil row sum

template<class MeshType>
inline scalar rowSum
(
    const meshDirection<stencil,MeshType>& st,
    const labelVector ijk
)
{
    scalar r(st(ijk).center());

    for (int c = 1; c < stencil::nComponents; c++)
        r += st(ijk)[c];

    return r;
}

template<class MeshType>
inline scalar rowSum
(
    const meshDirection<stencil,MeshType>& st,
    const label& i,
    const label& j,
    const label& k
)
{
    return rowSum(st,labelVector(i,j,k));
}

// Row multiplication

template<class Type, class MeshType>
inline Type rowProduct
(
    const meshDirection<stencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const label& i,
    const label& j,
    const label& k
)
{
    return
        st(i,j,k).center() * s(i,  j,  k  )
      + st(i,j,k).left()   * s(i-1,j,  k  )
      + st(i,j,k).right()  * s(i+1,j,  k  )
      + st(i,j,k).bottom() * s(i,  j-1,k  )
      + st(i,j,k).top()    * s(i,  j+1,k  )
      + st(i,j,k).aft()    * s(i,  j,  k-1)
      + st(i,j,k).fore()   * s(i,  j,  k+1);
}

template<class Type, class MeshType>
inline Type rowProduct
(
    const meshDirection<stencil,MeshType>& st,
    const Type& v,
    const label& i,
    const label& j,
    const label& k
)
{
    return rowSum(st,i,j,k)*v;
}

template<class Type, class MeshType>
inline Type rowProduct
(
    const meshDirection<stencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const labelVector& ijk
)
{
    return rowProduct(st(ijk),s,ijk.x(),ijk.y(),ijk.z());
}

template<class Type, class MeshType>
inline Type rowProduct
(
    const meshDirection<stencil,MeshType>& st,
    const Type& v,
    const labelVector& ijk
)
{
    return rowProduct(st(ijk),v,ijk.x(),ijk.y(),ijk.z());
}

#define STENCIL_PRODUCT(Type)                                                  \
template<>                                                                     \
class stencilProduct<stencil,Type>{ public: typedef Type type; };              \
template<>                                                                     \
class stencilProduct<Type,stencil>{ public: typedef Type type; };

STENCIL_PRODUCT(label)
STENCIL_PRODUCT(scalar)
STENCIL_PRODUCT(vector)
STENCIL_PRODUCT(tensor)
STENCIL_PRODUCT(symmTensor)
STENCIL_PRODUCT(sphericalTensor)
STENCIL_PRODUCT(diagTensor)

#undef STENCIL_PRODUCT

}

}

}

#endif
