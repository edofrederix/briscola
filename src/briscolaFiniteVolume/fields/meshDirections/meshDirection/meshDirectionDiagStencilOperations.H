#ifndef meshDirectionDiagStencilOperations_H
#define meshDirectionDiagStencilOperations_H

#include "diagStencil.H"
#include "meshDirectionStencilProductFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Stencil row sum

template<class MeshType>
inline scalar rowSum
(
    const meshDirection<diagStencil,MeshType>& st,
    const label& i,
    const label& j,
    const label& k
)
{
    return st(i,j,k).center();
}

template<class MeshType>
inline scalar rowSum
(
    const meshDirection<diagStencil,MeshType>& st,
    const labelVector ijk
)
{
    return st(ijk).center();
}

// Lower stencil row sum

template<class MeshType>
inline scalar lowerRowSum
(
    const meshDirection<diagStencil,MeshType>& st,
    const label& i,
    const label& j,
    const label& k
)
{
    return Zero;
}

template<class MeshType>
inline scalar lowerRowSum
(
    const meshDirection<diagStencil,MeshType>& st,
    const labelVector ijk
)
{
    return Zero;
}

// Upper stencil row sum

template<class MeshType>
inline scalar upperRowSum
(
    const meshDirection<diagStencil,MeshType>& st,
    const label& i,
    const label& j,
    const label& k
)
{
    return Zero;
}

template<class MeshType>
inline scalar upperRowSum
(
    const meshDirection<diagStencil,MeshType>& st,
    const labelVector ijk
)
{
    return Zero;
}

// Boundary source

template<class Type, class MeshType>
inline Type boundarySource
(
    const meshDirection<diagStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const label& faceNum,
    const labelVector& bo,
    const labelVector& ijk
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type boundarySource
(
    const meshDirection<diagStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const label& faceNum,
    const labelVector& bo,
    const label& i,
    const label& j,
    const label& k
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type boundarySource
(
    const meshDirection<diagStencil,MeshType>& st,
    const Type& v,
    const labelVector& bo,
    const label& faceNum,
    const labelVector& ijk
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type boundarySource
(
    const meshDirection<diagStencil,MeshType>& st,
    const Type& v,
    const label& faceNum,
    const labelVector& bo,
    const label& i,
    const label& j,
    const label& k
)
{
    return Zero;
}

// Row product

template<class Type, class MeshType>
inline Type rowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const label& i,
    const label& j,
    const label& k
)
{
    return st(i,j,k).center()*s(i,j,k);
}

template<class Type, class MeshType>
inline Type rowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const Type& v,
    const label& i,
    const label& j,
    const label& k
)
{
    return st(i,j,k).center()*v;
}

template<class Type, class MeshType>
inline Type rowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const labelVector& ijk
)
{
    return st(ijk)*s(ijk);
}

template<class Type, class MeshType>
inline Type rowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const Type& v,
    const labelVector& ijk
)
{
    return st(ijk)*v;
}

// Lower row product

template<class Type, class MeshType>
inline Type lowerRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const label& i,
    const label& j,
    const label& k
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type lowerRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const Type& v,
    const label& i,
    const label& j,
    const label& k
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type lowerRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const labelVector& ijk
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type lowerRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const Type& v,
    const labelVector& ijk
)
{
    return Zero;
}

// Upper row product

template<class Type, class MeshType>
inline Type upperRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const label& i,
    const label& j,
    const label& k
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type upperRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const Type& v,
    const label& i,
    const label& j,
    const label& k
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type upperRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const labelVector& ijk
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type upperRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const Type& v,
    const labelVector& ijk
)
{
    return Zero;
}

// Off-diagonal row product

template<class Type, class MeshType>
inline Type offDiagRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const label& i,
    const label& j,
    const label& k
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type offDiagRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const Type& v,
    const label& i,
    const label& j,
    const label& k
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type offDiagRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const meshDirection<Type,MeshType>& s,
    const labelVector& ijk
)
{
    return Zero;
}

template<class Type, class MeshType>
inline Type offDiagRowProduct
(
    const meshDirection<diagStencil,MeshType>& st,
    const Type& v,
    const labelVector& ijk
)
{
    return Zero;
}

#define DIAGSTENCIL_PRODUCT(Type)                                              \
template<>                                                                     \
class stencilProduct<diagStencil,Type>{ public: typedef Type type; };          \
template<>                                                                     \
class stencilProduct<Type,diagStencil>{ public: typedef Type type; };

DIAGSTENCIL_PRODUCT(label)
DIAGSTENCIL_PRODUCT(scalar)
DIAGSTENCIL_PRODUCT(vector)
DIAGSTENCIL_PRODUCT(tensor)
DIAGSTENCIL_PRODUCT(symmTensor)
DIAGSTENCIL_PRODUCT(sphericalTensor)
DIAGSTENCIL_PRODUCT(diagTensor)

#undef STENCIL_PRODUCT

}

}

}

#endif
