#define TEMPLATE template<class Type, class MeshType>
#include "meshDirectionFunctionsM.H"
#include "UPstream.H"

#define SCALARPRODTYPE typename scalarProduct<Type,Type>::type

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void mag(meshDirection<SCALARPRODTYPE,MeshType>&, const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshDirection<SCALARPRODTYPE,MeshType>>
mag(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshDirection<SCALARPRODTYPE,MeshType>>
mag(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
void cmptMax
(
    meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptMax(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptMax(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
void cmptMin
(
    meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptMin(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptMin(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
void cmptAv
(
    meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptAv(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptAv(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
void
cmptMag(meshDirection<Type,MeshType>&, const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshDirection<Type,MeshType>>
cmptMag(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshDirection<Type,MeshType>>
cmptMag(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
Type max(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
Type max(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
Type min(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
Type min(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
Type sum(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
Type sum(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
Type average(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
Type average(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
SCALARPRODTYPE maxMagSqr(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
SCALARPRODTYPE maxMagSqr(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
SCALARPRODTYPE minMagSqr(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
SCALARPRODTYPE minMagSqr(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
SCALARPRODTYPE sumProd
(
    const meshDirection<Type,MeshType>&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class MeshType>
SCALARPRODTYPE sumProd
(
    const tmp<meshDirection<Type,MeshType>>&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class MeshType>
SCALARPRODTYPE sumProd
(
    const meshDirection<Type,MeshType>&,
    const tmp<meshDirection<Type,MeshType>>&
);

template<class Type, class MeshType>
SCALARPRODTYPE sumProd
(
    const tmp<meshDirection<Type,MeshType>>&,
    const tmp<meshDirection<Type,MeshType>>&
);

template<class Type, class MeshType>
Type sumCmptProd
(
    const meshDirection<Type,MeshType>&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class MeshType>
Type sumCmptProd
(
    const tmp<meshDirection<Type,MeshType>>&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class MeshType>
Type sumCmptProd
(
    const meshDirection<Type,MeshType>&,
    const tmp<meshDirection<Type,MeshType>>&
);

template<class Type, class MeshType>
Type sumCmptProd
(
    const tmp<meshDirection<Type,MeshType>>&,
    const tmp<meshDirection<Type,MeshType>>&
);

template<class Type, class MeshType>
SCALARPRODTYPE sumMag(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
SCALARPRODTYPE sumMag(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
Type sumCmptMag(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
Type sumCmptMag(const tmp<meshDirection<Type,MeshType>>&);

#define G_UNARY_FUNCTION(ReturnType, gFunc)                                     \
                                                                                \
template<class Type, class MeshType>                                            \
ReturnType gFunc                                                                \
(                                                                               \
    const meshDirection<Type,MeshType>&,                                        \
    const label comm = UPstream::worldComm                                      \
);                                                                              \
                                                                                \
template<class Type, class MeshType>                                            \
ReturnType gFunc                                                                \
(                                                                               \
    const tmp<meshDirection<Type,MeshType>>&,                                   \
    const label comm = UPstream::worldComm                                      \
);

G_UNARY_FUNCTION(Type, gMax)
G_UNARY_FUNCTION(Type, gMin)
G_UNARY_FUNCTION(Type, gSum)
G_UNARY_FUNCTION(Type, gAverage)
G_UNARY_FUNCTION(SCALARPRODTYPE, gMaxMagSqr)
G_UNARY_FUNCTION(SCALARPRODTYPE, gMinMagSqr)
G_UNARY_FUNCTION(SCALARPRODTYPE, gSumMag)
G_UNARY_FUNCTION(Type, gSumCmptMag)

#undef G_UNARY_FUNCTION

BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)

UNARY_OPERATOR(Type, Type, -, negate)

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                          \
                                                                                       \
template<class Type1, class Type2, class MeshType>                                     \
void OpFunc                                                                            \
(                                                                                      \
    meshDirection<typename product<Type1, Type2>::type,MeshType>&,                     \
    const meshDirection<Type1,MeshType>&,                                              \
    const meshDirection<Type2,MeshType>&                                               \
);                                                                                     \
                                                                                       \
template<class Type1, class Type2, class MeshType>                                     \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>                      \
operator Op                                                                            \
(                                                                                      \
    const meshDirection<Type1,MeshType>&,                                              \
    const meshDirection<Type2,MeshType>&                                               \
);                                                                                     \
                                                                                       \
template<class Type1, class Type2, class MeshType>                                     \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>                      \
operator Op                                                                            \
(                                                                                      \
    const meshDirection<Type1,MeshType>&,                                              \
    const tmp<meshDirection<Type2,MeshType>>&                                          \
);                                                                                     \
                                                                                       \
template<class Type1, class Type2, class MeshType>                                     \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>                      \
operator Op                                                                            \
(                                                                                      \
    const tmp<meshDirection<Type1,MeshType>>&,                                         \
    const meshDirection<Type2,MeshType>&                                               \
);                                                                                     \
                                                                                       \
template<class Type1, class Type2, class MeshType>                                     \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>                      \
operator Op                                                                            \
(                                                                                      \
    const tmp<meshDirection<Type1,MeshType>>&,                                         \
    const tmp<meshDirection<Type2,MeshType>>&                                          \
);                                                                                     \
                                                                                       \
/* VectorSpace */                                                                      \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
void OpFunc                                                                            \
(                                                                                      \
    meshDirection<typename product<Type, Form>::type,MeshType>&,                       \
    const meshDirection<Type,MeshType>&,                                               \
    const VectorSpace<Form,Cmpt,nCmpt>&                                                \
);                                                                                     \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                        \
operator Op                                                                            \
(                                                                                      \
    const meshDirection<Type,MeshType>&,                                               \
    const VectorSpace<Form,Cmpt,nCmpt>&                                                \
);                                                                                     \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                        \
operator Op                                                                            \
(                                                                                      \
    const tmp<meshDirection<Type,MeshType>>&,                                          \
    const VectorSpace<Form,Cmpt,nCmpt>&                                                \
);                                                                                     \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
void OpFunc                                                                            \
(                                                                                      \
    meshDirection<typename product<Form, Type>::type,MeshType>&,                       \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                               \
    const meshDirection<Type,MeshType>&                                                \
);                                                                                     \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                        \
operator Op                                                                            \
(                                                                                      \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                               \
    const meshDirection<Type,MeshType>&                                                \
);                                                                                     \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                        \
operator Op                                                                            \
(                                                                                      \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                               \
    const tmp<meshDirection<Type,MeshType>>&                                           \
);                                                                                     \
                                                                                       \
/* CellSpace */                                                                        \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
void OpFunc                                                                            \
(                                                                                      \
    meshDirection<typename product<Type, Form>::type,MeshType>&,                       \
    const meshDirection<Type,MeshType>&,                                               \
    const CellSpace<Form,Cmpt,nCmpt>&                                                  \
);                                                                                     \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                        \
operator Op                                                                            \
(                                                                                      \
    const meshDirection<Type,MeshType>&,                                               \
    const CellSpace<Form,Cmpt,nCmpt>&                                                  \
);                                                                                     \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                        \
operator Op                                                                            \
(                                                                                      \
    const tmp<meshDirection<Type,MeshType>>&,                                          \
    const CellSpace<Form,Cmpt,nCmpt>&                                                  \
);                                                                                     \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
void OpFunc                                                                            \
(                                                                                      \
    meshDirection<typename product<Form, Type>::type,MeshType>&,                       \
    const CellSpace<Form,Cmpt,nCmpt>&,                                                 \
    const meshDirection<Type,MeshType>&                                                \
);                                                                                     \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                        \
operator Op                                                                            \
(                                                                                      \
    const CellSpace<Form,Cmpt,nCmpt>&,                                                 \
    const meshDirection<Type,MeshType>&                                                \
);                                                                                     \
                                                                                       \
template                                                                               \
<                                                                                      \
    class Type,                                                                        \
    class Form,                                                                        \
    class Cmpt,                                                                        \
    direction nCmpt,                                                                   \
    class MeshType                                                                     \
>                                                                                      \
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                        \
operator Op                                                                            \
(                                                                                      \
    const CellSpace<Form,Cmpt,nCmpt>&,                                                 \
    const tmp<meshDirection<Type,MeshType>>&                                           \
);

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)
PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR

}

}

}

#undef SCALARPRODTYPE

#include "undefBlockFunctionsM.H"

#include "scalarMeshDirection.H"
#include "faceScalarMeshDirection.H"
#include "faceVectorMeshDirection.H"
#include "diagStencilMeshDirection.H"
#include "stencilMeshDirection.H"
