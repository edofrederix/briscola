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

// Full stencil conversion

template<class MeshType>
inline stencil fullStencil
(
    const meshDirection<diagStencil,MeshType>& st,
    const label& i,
    const label& j,
    const label& k
)
{
    return static_cast<const stencil&>(st(i,j,k));
}

template<class MeshType>
inline stencil fullStencil
(
    const meshDirection<diagStencil,MeshType>& st,
    const labelVector ijk
)
{
    return static_cast<const stencil&>(st(ijk));
}

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

// Matrix row multiplication

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
