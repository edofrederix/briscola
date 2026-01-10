#define TEMPLATE template<class Type, class MeshType>
#include "meshLevelFunctionsM.H"
#include "UPstream.H"

// Scalar return type must be deduced because of cell space
#define SCALARPRODTYPE typename scalarProduct<Type,Type>::type

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void mag(meshLevel<SCALARPRODTYPE,MeshType>&, const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<SCALARPRODTYPE,MeshType>> mag(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<SCALARPRODTYPE,MeshType>> mag(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
void cmptMax
(
    meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>&,
    const meshLevel<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptMax(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptMax(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
void cmptMin
(
    meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>&,
    const meshLevel<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptMin(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptMin(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
void cmptAv
(
    meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>&,
    const meshLevel<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptAv(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptAv(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
void cmptMag(meshLevel<Type,MeshType>&, const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptMag(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptMag(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
void cmptSqr(meshLevel<Type,MeshType>&, const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptSqr(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptSqr(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
void cmptSqrt(meshLevel<Type,MeshType>&, const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptSqrt(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptSqrt(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
List<Type> max(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
List<Type> max(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
List<Type> min(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
List<Type> min(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
List<Type> sum(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
List<Type> sum(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
List<Type> average(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
List<Type> average(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
List<SCALARPRODTYPE> maxMagSqr(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
List<SCALARPRODTYPE> maxMagSqr(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
List<SCALARPRODTYPE> minMagSqr(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
List<SCALARPRODTYPE> minMagSqr(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
List<SCALARPRODTYPE> sumProd
(
    const meshLevel<Type,MeshType>&,
    const meshLevel<Type,MeshType>&
);

template<class Type, class MeshType>
List<SCALARPRODTYPE> sumProd
(
    const tmp<meshLevel<Type,MeshType>>&,
    const meshLevel<Type,MeshType>&
);

template<class Type, class MeshType>
List<SCALARPRODTYPE> sumProd
(
    const meshLevel<Type,MeshType>&,
    const tmp<meshLevel<Type,MeshType>>&
);

template<class Type, class MeshType>
List<SCALARPRODTYPE> sumProd
(
    const tmp<meshLevel<Type,MeshType>>&,
    const tmp<meshLevel<Type,MeshType>>&
);

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const meshLevel<Type,MeshType>&,
    const meshLevel<Type,MeshType>&
);

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const tmp<meshLevel<Type,MeshType>>&,
    const meshLevel<Type,MeshType>&
);

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const meshLevel<Type,MeshType>&,
    const tmp<meshLevel<Type,MeshType>>&
);

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const tmp<meshLevel<Type,MeshType>>&,
    const tmp<meshLevel<Type,MeshType>>&
);

template<class Type, class MeshType>
List<SCALARPRODTYPE> sumMag(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
List<SCALARPRODTYPE> sumMag(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
List<Type> sumCmptMag(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
List<Type> sumCmptMag(const tmp<meshLevel<Type,MeshType>>&);

#define G_UNARY_FUNCTION(ReturnType, gFunc)                                    \
                                                                               \
template<class Type, class MeshType>                                           \
List<ReturnType> gFunc(const meshLevel<Type,MeshType>&);                       \
                                                                               \
template<class Type, class MeshType>                                           \
List<ReturnType> gFunc(const tmp<meshLevel<Type,MeshType>>&);

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
    meshLevel<typename product<Type1, Type2>::type,MeshType>&,                         \
    const meshLevel<Type1,MeshType>&,                                                  \
    const meshLevel<Type2,MeshType>&                                                   \
);                                                                                     \
                                                                                       \
template<class Type1, class Type2, class MeshType>                                     \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                          \
operator Op                                                                            \
(                                                                                      \
    const meshLevel<Type1,MeshType>&,                                                  \
    const meshLevel<Type2,MeshType>&                                                   \
);                                                                                     \
                                                                                       \
template<class Type1, class Type2, class MeshType>                                     \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                          \
operator Op                                                                            \
(                                                                                      \
    const meshLevel<Type1,MeshType>&,                                                  \
    const tmp<meshLevel<Type2,MeshType>>&                                              \
);                                                                                     \
                                                                                       \
template<class Type1, class Type2, class MeshType>                                     \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                          \
operator Op                                                                            \
(                                                                                      \
    const tmp<meshLevel<Type1,MeshType>>&,                                             \
    const meshLevel<Type2,MeshType>&                                                   \
);                                                                                     \
                                                                                       \
template<class Type1, class Type2, class MeshType>                                     \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                          \
operator Op                                                                            \
(                                                                                      \
    const tmp<meshLevel<Type1,MeshType>>&,                                             \
    const tmp<meshLevel<Type2,MeshType>>&                                              \
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
    meshLevel<typename product<Type, Form>::type,MeshType>&,                           \
    const meshLevel<Type,MeshType>&,                                                   \
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
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const meshLevel<Type,MeshType>&,                                                   \
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
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const tmp<meshLevel<Type,MeshType>>&,                                              \
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
    meshLevel<typename product<Form, Type>::type,MeshType>&,                           \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                               \
    const meshLevel<Type,MeshType>&                                                    \
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
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                               \
    const meshLevel<Type,MeshType>&                                                    \
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
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                               \
    const tmp<meshLevel<Type,MeshType>>&                                               \
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
    meshLevel<typename product<Type, Form>::type,MeshType>&,                           \
    const meshLevel<Type,MeshType>&,                                                   \
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
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const meshLevel<Type,MeshType>&,                                                   \
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
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const tmp<meshLevel<Type,MeshType>>&,                                              \
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
    meshLevel<typename product<Form, Type>::type,MeshType>&,                           \
    const CellSpace<Form,Cmpt,nCmpt>&,                                                 \
    const meshLevel<Type,MeshType>&                                                    \
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
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const CellSpace<Form,Cmpt,nCmpt>&,                                                 \
    const meshLevel<Type,MeshType>&                                                    \
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
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const CellSpace<Form,Cmpt,nCmpt>&,                                                 \
    const tmp<meshLevel<Type,MeshType>>&                                               \
);                                                                                     \
                                                                                       \
/* Lists */                                                                            \
                                                                                       \
template<class Type, class Form, class MeshType>                                       \
void OpFunc                                                                            \
(                                                                                      \
    meshLevel<typename product<Type, Form>::type,MeshType>&,                           \
    const meshLevel<Type,MeshType>&,                                                   \
    const List<Form>&                                                                  \
);                                                                                     \
                                                                                       \
template<class Type, class Form, class MeshType>                                       \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const meshLevel<Type,MeshType>&,                                                   \
    const List<Form>&                                                                  \
);                                                                                     \
                                                                                       \
template<class Type, class Form, class MeshType>                                       \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const tmp<meshLevel<Type,MeshType>>&,                                              \
    const List<Form>&                                                                  \
);                                                                                     \
                                                                                       \
template<class Form, class Type, class MeshType>                                       \
void OpFunc                                                                            \
(                                                                                      \
    meshLevel<typename product<Form, Type>::type,MeshType>&,                           \
    const List<Form>&,                                                                 \
    const meshLevel<Type,MeshType>&                                                    \
);                                                                                     \
                                                                                       \
template<class Form, class Type, class MeshType>                                       \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const List<Form>&,                                                                 \
    const meshLevel<Type,MeshType>&                                                    \
);                                                                                     \
                                                                                       \
template<class Form, class Type, class MeshType>                                       \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                            \
operator Op                                                                            \
(                                                                                      \
    const List<Form>&,                                                                 \
    const tmp<meshLevel<Type,MeshType>>&                                               \
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

#include "scalarMeshLevel.H"
#include "labelMeshLevel.H"
#include "vectorMeshLevel.H"
#include "tensorMeshLevel.H"
#include "symmTensorMeshLevel.H"
#include "sphericalTensorMeshLevel.H"
#include "diagTensorMeshLevel.H"
#include "faceScalarMeshLevel.H"
#include "edgeScalarMeshLevel.H"
#include "vertexScalarMeshLevel.H"
#include "faceVectorMeshLevel.H"
#include "edgeVectorMeshLevel.H"
#include "vertexVectorMeshLevel.H"
#include "diagStencilMeshLevel.H"
#include "stencilMeshLevel.H"
