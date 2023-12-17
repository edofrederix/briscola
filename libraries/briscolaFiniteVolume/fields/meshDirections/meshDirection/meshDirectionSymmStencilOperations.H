#ifndef meshDirectionSymmStencilOperations_H
#define meshDirectionSymmStencilOperations_H

#include "symmStencil.H"
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
    const meshDirection<symmStencil,MeshType>& st,
    const labelVector ijk
)
{
    return stencil
    (
        st(ijk      ).center(),
        st(ijk      ).left(),
        st(ijk+unitX).left(),
        st(ijk      ).bottom(),
        st(ijk+unitY).bottom(),
        st(ijk      ).aft(),
        st(ijk+unitZ).aft()
    );
}

template<class MeshType>
inline stencil fullStencil
(
    const meshDirection<symmStencil,MeshType>& st,
    const label& i,
    const label& j,
    const label& k
)
{
    return fullStencil(st,labelVector(i,j,k));
}

// Stencil row sum

template<class MeshType>
inline scalar rowSum
(
    const meshDirection<symmStencil,MeshType>& st,
    const labelVector ijk
)
{
    scalar r(0);

    for (int c = 0; c < symmStencil::nComponents; c++)
        r += st(ijk)[c];

    r += st(ijk+unitX).left();
    r += st(ijk+unitY).bottom();
    r += st(ijk+unitZ).aft();

    return r;
}

template<class MeshType>
inline scalar rowSum
(
    const meshDirection<symmStencil,MeshType>& st,
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
inline Type rowProduct
(
    const meshDirection<symmStencil,MeshType>& st,
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
    const meshDirection<symmStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const labelVector& ijk
)
{
    return rowProduct(st,s,ijk.x(),ijk.y(),ijk.z());
}

template<class Type, class MeshType>
inline Type rowProduct
(
    const meshDirection<symmStencil,MeshType>& st,
    const Type& v,
    const labelVector& ijk
)
{
    return rowProduct(st,v,ijk.x(),ijk.y(),ijk.z());
}

#define STENCIL_PRODUCT(Type)                                                  \
template<>                                                                     \
class stencilProduct<symmStencil,Type>{ public: typedef Type type; };          \
template<>                                                                     \
class stencilProduct<Type,symmStencil>{ public: typedef Type type; };

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
