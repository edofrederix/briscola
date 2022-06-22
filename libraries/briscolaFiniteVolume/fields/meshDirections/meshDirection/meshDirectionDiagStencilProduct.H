#ifndef meshDirectionDiagStencilProduct_H
#define meshDirectionDiagStencilProduct_H

#include "diagStencil.H"
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
inline Type Amul
(
    const meshDirection<Type,MeshType>& s,
    const meshDirection<diagStencil,MeshType>& st,
    const label& i,
    const label& j,
    const label& k
)
{
    return st(i,j,k).center()*s(i,j,k);
}

template<class Type, class MeshType>
inline Type Amul
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

#define DIAGSTENCIL_PRODUCT(Type)                                               \
template<>                                                                      \
class stencilProduct<diagStencil,Type>{ public: typedef Type type; };           \
template<>                                                                      \
class stencilProduct<Type,diagStencil>{ public: typedef Type type; };

DIAGSTENCIL_PRODUCT(label)
DIAGSTENCIL_PRODUCT(scalar)
DIAGSTENCIL_PRODUCT(vector)
DIAGSTENCIL_PRODUCT(tensor)
DIAGSTENCIL_PRODUCT(symmTensor)
DIAGSTENCIL_PRODUCT(sphericalTensor)
DIAGSTENCIL_PRODUCT(diagTensor)

}

}

}

#endif
