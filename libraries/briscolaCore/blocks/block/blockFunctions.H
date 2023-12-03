#define TEMPLATE template<class Type>
#include "blockFunctionsM.H"
#include "UPstream.H"

// Scalar return type must be deduced because of cell space
#define SCALARPRODTYPE typename scalarProduct<Type,Type>::type

namespace Foam
{

namespace briscola
{

template<class Type>
void mag(block<SCALARPRODTYPE>&, const block<Type>&);

template<class Type>
tmp<block<SCALARPRODTYPE>> mag(const block<Type>&);

template<class Type>
tmp<block<SCALARPRODTYPE>> mag(const tmp<block<Type>>&);

template<class Type>
void cmptMax(block<typename block<Type>::cmptType>&, const block<Type>&);

template<class Type>
tmp<block<typename block<Type>::cmptType>>
cmptMax(const block<Type>&);

template<class Type>
tmp<block<typename block<Type>::cmptType>>
cmptMax(const tmp<block<Type>>&);

template<class Type>
void cmptMin(block<typename block<Type>::cmptType>&, const block<Type>&);

template<class Type>
tmp<block<typename block<Type>::cmptType>>
cmptMin(const block<Type>&);

template<class Type>
tmp<block<typename block<Type>::cmptType>>
cmptMin(const tmp<block<Type>>&);

template<class Type>
void cmptAv(block<typename block<Type>::cmptType>&, const block<Type>&);

template<class Type>
tmp<block<typename block<Type>::cmptType>>
cmptAv(const block<Type>&);

template<class Type>
tmp<block<typename block<Type>::cmptType>>
cmptAv(const tmp<block<Type>>&);

template<class Type>
void cmptMag(block<Type>&, const block<Type>&);

template<class Type>
tmp<block<Type>> cmptMag(const block<Type>&);

template<class Type>
tmp<block<Type>> cmptMag(const tmp<block<Type>>&);


template<class Type>
Type max(const block<Type>&);

template<class Type>
Type max(const tmp<block<Type>>&);

template<class Type>
Type min(const block<Type>&);

template<class Type>
Type min(const tmp<block<Type>>&);

template<class Type>
Type sum(const block<Type>&);

template<class Type>
Type sum(const tmp<block<Type>>&);

template<class Type>
Type average(const block<Type>&);

template<class Type>
Type average(const tmp<block<Type>>&);


template<class Type>
SCALARPRODTYPE maxMagSqr(const block<Type>&);

template<class Type>
SCALARPRODTYPE maxMagSqr(const tmp<block<Type>>&);

template<class Type>
SCALARPRODTYPE minMagSqr(const block<Type>&);

template<class Type>
SCALARPRODTYPE minMagSqr(const tmp<block<Type>>&);

template<class Type>
SCALARPRODTYPE sumProd(const block<Type>&, const block<Type>&);

template<class Type>
SCALARPRODTYPE sumProd(const tmp<block<Type>>&, const block<Type>&);

template<class Type>
SCALARPRODTYPE sumProd(const block<Type>&, const tmp<block<Type>>&);

template<class Type>
SCALARPRODTYPE sumProd(const tmp<block<Type>>&, const tmp<block<Type>>&);


template<class Type>
Type sumCmptProd(const block<Type>&, const block<Type>&);

template<class Type>
Type sumCmptProd(const tmp<block<Type>>&, const block<Type>&);

template<class Type>
Type sumCmptProd(const block<Type>&, const tmp<block<Type>>&);

template<class Type>
Type sumCmptProd(const tmp<block<Type>>&, const tmp<block<Type>>&);

template<class Type>
typename powProduct<Type,2>::type sumSqr(const block<Type>&);

template<class Type>
typename powProduct<Type,2>::type sumSqr(const tmp<block<Type>>&);

template<class Type>
SCALARPRODTYPE sumMag(const block<Type>&);

template<class Type>
SCALARPRODTYPE sumMag(const tmp<block<Type>>&);

template<class Type>
Type sumCmptMag(const block<Type>&);

template<class Type>
Type sumCmptMag(const tmp<block<Type>>&);

#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                        \
                                                                                \
template<class Type>                                                            \
ReturnType gFunc                                                                \
(                                                                               \
    const block<Type>&,                                                         \
    const label comm = UPstream::worldComm                                      \
);                                                                              \
                                                                                \
template<class Type>                                                            \
ReturnType gFunc                                                                \
(                                                                               \
    const tmp<block<Type>>&,                                                    \
    const label comm = UPstream::worldComm                                      \
);

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(Type, gAverage, average, average)
G_UNARY_FUNCTION(SCALARPRODTYPE, gMaxMagSqr, maxMagSqr, maxMagSqr)
G_UNARY_FUNCTION(SCALARPRODTYPE, gMinMagSqr, minMagSqr, minMagSqr)
G_UNARY_FUNCTION(SCALARPRODTYPE, gSumMag, sumMag, sum)
G_UNARY_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)

#undef G_UNARY_FUNCTION

template<class Type>
typename powProduct<Type,2>::type gSumSqr
(
    const block<Type>&,
    const label comm = UPstream::worldComm
);

template<class Type>
typename powProduct<Type,2>::type gSumSqr
(
    const tmp<block<Type>>&,
    const label comm = UPstream::worldComm
);

// template<Type>
//
//     block<arg1> = arg4(block<arg2>, block<arg3>)

BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

// template<Type>
//
//     block<arg1> = arg4(arg2, block<arg3>)
//     block<arg1> = arg4(block<arg2>, arg3)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)

// template<Type>
//
//     block<arg1> = arg3 block<arg2>

UNARY_OPERATOR(Type, Type, -, negate)

// template<Type>
//
//     block<arg1> = block<arg2> arg4 block<arg3>

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

// template<Type>
//
//     block<arg1> = arg2 arg4 block<arg3>

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)

// template<Type>
//
//     block<arg1> = block<arg2> arg4 arg3

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type1, Type2>::type>&,                              \
    const block<Type1>&,                                                       \
    const block<Type2>&                                                        \
);                                                                             \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const block<Type1>&, const block<Type2>&);                         \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const block<Type1>&, const tmp<block<Type2>>&);                    \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<block<Type1>>&, const block<Type2>&);                    \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<block<Type1>>&, const tmp<block<Type2>>&);               \
                                                                               \
/* VectorSpace */                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type, Form>::type>&,                                \
    const block<Type>&,                                                        \
    const VectorSpace<Form,Cmpt,nCmpt>&                                        \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op(const block<Type>&, const VectorSpace<Form,Cmpt,nCmpt>&);          \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op(const tmp<block<Type>>&, const VectorSpace<Form,Cmpt,nCmpt>&);     \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Form, Type>::type>&,                                \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                       \
    const block<Type>&                                                         \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>&, const block<Type>&);          \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>&,const tmp<block<Type>>&);      \
                                                                               \
/* CellSpace */                                                                \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type, Form>::type>&,                                \
    const block<Type>&,                                                        \
    const CellSpace<Form,Cmpt,nCmpt>&                                          \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op(const block<Type>&, const CellSpace<Form,Cmpt,nCmpt>&);            \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op(const tmp<block<Type>>&, const CellSpace<Form,Cmpt,nCmpt>&);       \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Form, Type>::type>&,                                \
    const CellSpace<Form,Cmpt,nCmpt>&,                                         \
    const block<Type>&                                                         \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op(const CellSpace<Form,Cmpt,nCmpt>&, const block<Type>&);            \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op(const CellSpace<Form,Cmpt,nCmpt>&,const tmp<block<Type>>&);

// template<Type1, Type2>, VS = VectorSpace, CS = CellSpace
//
//     block<arg1<Type1,Type2>> = block<Type1> arg2 block<Type2>
//     block<arg1<Type1,Type2>> = VS<Type1> arg2 block<Type2>
//     block<arg1<Type1,Type2>> = block<Type1> arg2 VS<Type2>
//     block<arg1<Type1,Type2>> = CS<Type1> arg2 block<Type2>
//     block<arg1<Type1,Type2>> = block<Type1> arg2 CS<Type2>
//
// Note: this does not define
//
//     block<arg1<Type1,Type2>> = Type1 arg2 block<Type2>
//     block<arg1<Type1,Type2>> = block<Type1> arg2 Type2
//
// because this generates unresolvable overloads.

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)
PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR

}

}

#undef SCALARPRODTYPE

#include "undefBlockFunctionsM.H"
#include "scalarBlock.H"
#include "faceScalarBlock.H"
#include "lowerFaceScalarBlock.H"
#include "edgeScalarBlock.H"
#include "vertexScalarBlock.H"
#include "faceVectorBlock.H"
#include "lowerFaceVectorBlock.H"
#include "edgeVectorBlock.H"
#include "vertexVectorBlock.H"
#include "diagStencilBlock.H"
#include "stencilBlock.H"
