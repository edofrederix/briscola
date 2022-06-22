#ifndef meshDirectionStencilProduct_H
#define meshDirectionStencilProduct_H

#include "stencil.H"
#include "meshDirectionStencilProductFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
inline Type Amul
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
inline Type Amul
(
    const meshDirection<Type,MeshType>& s,
    const meshDirection<stencil,MeshType>& st,
    const label& i,
    const label& j,
    const label& k
)
{
    return
        st(i,j,k).center() * s(i,  j,  k  )
      + st(i,j,k).right()  * s(i-1,j,  k  )
      + st(i,j,k).left()   * s(i+1,j,  k  )
      + st(i,j,k).top()    * s(i,  j-1,k  )
      + st(i,j,k).bottom() * s(i,  j+1,k  )
      + st(i,j,k).fore()   * s(i,  j,  k-1)
      + st(i,j,k).aft()    * s(i,  j,  k+1);
}

template<class Type, class MeshType>
inline Type Amul
(
    const meshDirection<stencil,MeshType>& st,
    const Type& v,
    const label& i,
    const label& j,
    const label& k
)
{
    return stencilSum(st(i,j,k))*v;
}

#define STENCIL_PRODUCT(Type)                                                   \
template<>                                                                      \
class stencilProduct<stencil,Type>{ public: typedef Type type; };               \
template<>                                                                      \
class stencilProduct<Type,stencil>{ public: typedef Type type; };

STENCIL_PRODUCT(label)
STENCIL_PRODUCT(scalar)
STENCIL_PRODUCT(vector)
STENCIL_PRODUCT(tensor)
STENCIL_PRODUCT(symmTensor)
STENCIL_PRODUCT(sphericalTensor)
STENCIL_PRODUCT(diagTensor)

}

}

}

#endif
