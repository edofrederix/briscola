#define TEMPLATE template<class Type, class MeshType>

#include "meshFieldFunctionsM.H"
#include "UPstream.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void mag
(
    meshField<scalar,MeshType>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
mag(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
mag(const tmp<meshField<Type,MeshType>>&);


template<class Type, class MeshType>
void cmptMax
(
    meshField<typename meshField<Type,MeshType>::cmptType,MeshType>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptMax
(
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptMax
(
    const tmp<meshField<Type,MeshType>>&
);

template<class Type, class MeshType>
void cmptMin
(
    meshField<typename meshField<Type,MeshType>::cmptType,MeshType>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptMin
(
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptMin
(
    const tmp<meshField<Type,MeshType>>&
);

template<class Type, class MeshType>
void cmptAv
(
    meshField<typename meshField<Type,MeshType>::cmptType,MeshType>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptAv
(
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptAv
(
    const tmp<meshField<Type,MeshType>>&
);

template<class Type, class MeshType>
void cmptMag
(
    meshField<Type,MeshType>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
cmptMag
(
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
cmptMag
(
    const tmp<meshField<Type,MeshType>>&
);

template<class Type, class MeshType>
List<Type>
max(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
List<Type>
max(const tmp<meshField<Type,MeshType>>&);

template<class Type, class MeshType>
List<Type>
min(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
List<Type>
min(const tmp<meshField<Type,MeshType>>&);

template<class Type, class MeshType>
List<Type>
sum(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
List<Type>
sum(const tmp<meshField<Type,MeshType>>&);

template<class Type, class MeshType>
List<Type>
average(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
List<Type>
average(const tmp<meshField<Type,MeshType>>&);

template<class Type, class MeshType>
List<scalar>
maxMagSqr(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
List<scalar>
maxMagSqr(const tmp<meshField<Type,MeshType>>&);

template<class Type, class MeshType>
List<scalar>
minMagSqr(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
List<scalar>
minMagSqr(const tmp<meshField<Type,MeshType>>&);

template<class Type, class MeshType>
List<scalar> sumProd
(
    const meshField<Type,MeshType>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
List<scalar> sumProd
(
    const tmp<meshField<Type,MeshType>>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
List<scalar> sumProd
(
    const meshField<Type,MeshType>&,
    const tmp<meshField<Type,MeshType>>&
);

template<class Type, class MeshType>
List<scalar> sumProd
(
    const tmp<meshField<Type,MeshType>>&,
    const tmp<meshField<Type,MeshType>>&
);

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const meshField<Type,MeshType>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const tmp<meshField<Type,MeshType>>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const meshField<Type,MeshType>&,
    const tmp<meshField<Type,MeshType>>&
);

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const tmp<meshField<Type,MeshType>>&,
    const tmp<meshField<Type,MeshType>>&
);

template<class Type, class MeshType>
List<scalar>
sumMag(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
List<scalar>
sumMag(const tmp<meshField<Type,MeshType>>&);


template<class Type, class MeshType>
List<Type>
sumCmptMag(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
List<Type>
sumCmptMag(const tmp<meshField<Type,MeshType>>&);

#define G_UNARY_FUNCTION(ReturnType, gFunc)                                     \
                                                                                \
template<class Type, class MeshType>                                            \
List<ReturnType> gFunc                                                          \
(                                                                               \
    const meshField<Type,MeshType>&,                                            \
    const label comm = UPstream::worldComm                                      \
);                                                                              \
                                                                                \
template<class Type, class MeshType>                                            \
List<ReturnType> gFunc                                                          \
(                                                                               \
    const tmp<meshField<Type,MeshType>>&,                                       \
    const label comm = UPstream::worldComm                                      \
);                                                                              \

G_UNARY_FUNCTION(Type, gMax)
G_UNARY_FUNCTION(Type, gMin)
G_UNARY_FUNCTION(Type, gSum)
G_UNARY_FUNCTION(Type, gAverage)
G_UNARY_FUNCTION(scalar, gMaxMagSqr)
G_UNARY_FUNCTION(scalar, gMinMagSqr)
G_UNARY_FUNCTION(scalar, gSumMag)
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

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                   \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
void OpFunc                                                                     \
(                                                                               \
    meshField<typename product<Type1, Type2>::type,MeshType>&,                  \
    const meshField<Type1,MeshType>&,                                           \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshField<typename product<Type1, Type2>::type,MeshType>>                   \
operator Op                                                                     \
(                                                                               \
    const meshField<Type1,MeshType>&,                                           \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshField<typename product<Type1, Type2>::type,MeshType>>                   \
operator Op                                                                     \
(                                                                               \
    const meshField<Type1,MeshType>&,                                           \
    const tmp<meshField<Type2,MeshType>>&                                       \
);                                                                              \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshField<typename product<Type1, Type2>::type,MeshType>>                   \
operator Op                                                                     \
(                                                                               \
    const tmp<meshField<Type1,MeshType>>&,                                      \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshField<typename product<Type1, Type2>::type,MeshType>>                   \
operator Op                                                                     \
(                                                                               \
    const tmp<meshField<Type1,MeshType>>&,                                      \
    const tmp<meshField<Type2,MeshType>>&                                       \
);                                                                              \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
void OpFunc                                                                     \
(                                                                               \
    meshField<typename product<Type, Form>::type,MeshType>&,                    \
    const meshField<Type,MeshType>&,                                            \
    const VectorSpace<Form,Cmpt,nCmpt>&                                         \
);                                                                              \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshField<typename product<Type, Form>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const meshField<Type,MeshType>&,                                            \
    const VectorSpace<Form,Cmpt,nCmpt>&                                         \
);                                                                              \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshField<typename product<Type, Form>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const tmp<meshField<Type,MeshType>>&,                                       \
    const VectorSpace<Form,Cmpt,nCmpt>&                                         \
);                                                                              \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
void OpFunc                                                                     \
(                                                                               \
    meshField<typename product<Form,Type>::type,MeshType>&,                     \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                        \
    const meshField<Type,MeshType>&                                             \
);                                                                              \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshField<typename product<Form,Type>::type,MeshType>>                      \
operator Op                                                                     \
(                                                                               \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                        \
    const meshField<Type,MeshType>&                                             \
);                                                                              \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshField<typename product<Form,Type>::type,MeshType>>                      \
operator Op                                                                     \
(                                                                               \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                        \
    const tmp<meshField<Type,MeshType>>&                                        \
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

#include "undefBlockFunctionsM.H"
