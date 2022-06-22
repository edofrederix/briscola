#ifndef meshDirectionSymmStencilProduct_H
#define meshDirectionSsymmStencilProduct_H

#include "symmStencil.H"
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
    const meshDirection<symmStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const label& i,
    const label& j,
    const label& k
)
{
    return
        st(i,  j,  k  ).center() * s(i,  j,  k  )
      + st(i,  j,  k  ).left()   * s(i-1,j,  k  )
      + st(i+1,j,  k  ).left()   * s(i+1,j,  k  )
      + st(i,  j,  k  ).bottom() * s(i,  j-1,k  )
      + st(i,  j+1,k  ).bottom() * s(i,  j+1,k  )
      + st(i,  j,  k  ).aft()    * s(i,  j,  k-1)
      + st(i,  j,  k+1).aft()    * s(i,  j,  k+1);
}

template<class Type, class MeshType>
inline Type Amul
(
    const meshDirection<Type,MeshType>& s,
    const meshDirection<symmStencil,MeshType>& st,
    const label& i,
    const label& j,
    const label& k
)
{
    return
        st(i,  j,  k  ).center() * s(i,  j,  k  )
      + st(i+1,j,  k  ).left()   * s(i-1,j,  k  )
      + st(i,  j,  k  ).left()   * s(i+1,j,  k  )
      + st(i,  j+1,k  ).bottom() * s(i,  j-1,k  )
      + st(i,  j,  k  ).bottom() * s(i,  j+1,k  )
      + st(i,  j,  k+1).aft()    * s(i,  j,  k-1)
      + st(i,  j,  k  ).aft()    * s(i,  j,  k+1);
}

template<class Type, class MeshType>
inline Type Amul
(
    const meshDirection<symmStencil,MeshType>& st,
    const Type& v,
    const label& i,
    const label& j,
    const label& k
)
{
    return stencilSum(st(i,j,k))*v;
}

#define SYMMSTENCIL_PRODUCT(Type)                                               \
template<>                                                                      \
class stencilProduct<symmStencil,Type>{ public: typedef Type type; };           \
template<>                                                                      \
class stencilProduct<Type,symmStencil>{ public: typedef Type type; };

SYMMSTENCIL_PRODUCT(label)
SYMMSTENCIL_PRODUCT(scalar)
SYMMSTENCIL_PRODUCT(vector)
SYMMSTENCIL_PRODUCT(tensor)
SYMMSTENCIL_PRODUCT(symmTensor)
SYMMSTENCIL_PRODUCT(sphericalTensor)
SYMMSTENCIL_PRODUCT(diagTensor)

}

}

}

#endif
