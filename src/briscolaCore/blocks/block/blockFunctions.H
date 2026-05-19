#define TEMPLATE template<class Type, int P>
#include "blockFunctionsM.H"
#include "UPstream.H"

// Scalar return type must be deduced because of cell space
#define SCALARPRODTYPE typename scalarProduct<Type,Type>::type

namespace Foam
{

namespace briscola
{

template<class Type, int P>
void component
(
    block<typename block<Type,P>::cmptType, P>&,
    const block<Type,P>&,
    const direction
);

template<class Type, int P>
void T(block<Type,P>&, const block<Type,P>&);

template<class Type, int P>
void sqr
(
    block<typename outerProduct<Type, Type>::type, P>&,
    const block<Type,P>&
);

template<class Type, int P>
tmp<block<typename outerProduct<Type, Type>::type, P>>
sqr(const block<Type,P>&);

template<class Type, int P>
tmp<block<typename outerProduct<Type, Type>::type, P>>
sqr(const tmp<block<Type,P>>&);

template<class Type, int P>
void mag(block<SCALARPRODTYPE,P>&, const block<Type,P>&);

template<class Type, int P>
tmp<block<SCALARPRODTYPE,P>> mag(const block<Type,P>&);

template<class Type, int P>
tmp<block<SCALARPRODTYPE,P>> mag(const tmp<block<Type,P>>&);

template<class Type, int P>
void cmptMag(block<Type,P>&, const block<Type,P>&);

template<class Type, int P>
tmp<block<Type,P>> cmptMag(const block<Type,P>&);

template<class Type, int P>
tmp<block<Type,P>> cmptMag(const tmp<block<Type,P>>&);

template<class Type, int P>
void cmptSqr(block<Type,P>&, const block<Type,P>&);

template<class Type, int P>
tmp<block<Type,P>> cmptSqr(const block<Type,P>&);

template<class Type, int P>
tmp<block<Type,P>> cmptSqr(const tmp<block<Type,P>>&);

template<class Type, int P>
Type max(const block<Type,P>&);

template<class Type, int P>
Type max(const tmp<block<Type,P>>&);

template<class Type, int P>
Type min(const block<Type,P>&);

template<class Type, int P>
Type min(const tmp<block<Type,P>>&);

template<class Type, int P>
Type sum(const block<Type,P>&);

template<class Type, int P>
Type sum(const tmp<block<Type,P>>&);

template<class Type, int P>
Type average(const block<Type,P>&);

template<class Type, int P>
Type average(const tmp<block<Type,P>>&);

#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                       \
                                                                               \
template<class Type, int P>                                                    \
ReturnType gFunc                                                               \
(                                                                              \
    const block<Type,P>&,                                                      \
    const label comm = UPstream::worldComm                                     \
);                                                                             \
                                                                               \
template<class Type, int P>                                                    \
ReturnType gFunc                                                               \
(                                                                              \
    const tmp<block<Type,P>>&,                                                 \
    const label comm = UPstream::worldComm                                     \
);

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(Type, gAverage, average, average)

#undef G_UNARY_FUNCTION

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
template<class Type1, class Type2, int P>                                      \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type1, Type2>::type, P>&,                           \
    const block<Type1,P>&,                                                     \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
template<class Type1, class Type2, int P>                                      \
tmp<block<typename product<Type1, Type2>::type, P>>                            \
operator Op(const block<Type1,P>&, const block<Type2,P>&);                     \
                                                                               \
template<class Type1, class Type2, int P>                                      \
tmp<block<typename product<Type1, Type2>::type, P>>                            \
operator Op(const block<Type1,P>&, const tmp<block<Type2,P>>&);                \
                                                                               \
template<class Type1, class Type2, int P>                                      \
tmp<block<typename product<Type1, Type2>::type, P>>                            \
operator Op(const tmp<block<Type1,P>>&, const block<Type2,P>&);                \
                                                                               \
template<class Type1, class Type2, int P>                                      \
tmp<block<typename product<Type1, Type2>::type, P>>                            \
operator Op(const tmp<block<Type1,P>>&, const tmp<block<Type2,P>>&);           \
                                                                               \
/* VectorSpace */                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type, Form>::type, P>&,                             \
    const block<Type,P>&,                                                      \
    const VectorSpace<Form,Cmpt,nCmpt>&                                        \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Type, Form>::type, P>>                              \
operator Op(const block<Type,P>&, const VectorSpace<Form,Cmpt,nCmpt>&);        \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Type, Form>::type, P>>                              \
operator Op(const tmp<block<Type,P>>&, const VectorSpace<Form,Cmpt,nCmpt>&);   \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Form, Type>::type, P>&,                             \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                       \
    const block<Type,P>&                                                       \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Form, Type>::type, P>>                              \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>&, const block<Type,P>&);        \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Form, Type>::type, P>>                              \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>&, const tmp<block<Type,P>>&);   \
                                                                               \
/* CellSpace */                                                                \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type, Form>::type, P>&,                             \
    const block<Type,P>&,                                                      \
    const CellSpace<Form,Cmpt,nCmpt>&                                          \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Type, Form>::type, P>>                              \
operator Op(const block<Type,P>&, const CellSpace<Form,Cmpt,nCmpt>&);          \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Type, Form>::type, P>>                              \
operator Op(const tmp<block<Type,P>>&, const CellSpace<Form,Cmpt,nCmpt>&);     \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Form, Type>::type, P>&,                             \
    const CellSpace<Form,Cmpt,nCmpt>&,                                         \
    const block<Type,P>&                                                       \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Form, Type>::type, P>>                              \
operator Op(const CellSpace<Form,Cmpt,nCmpt>&, const block<Type,P>&);          \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Form, Type>::type, P>>                              \
operator Op(const CellSpace<Form,Cmpt,nCmpt>&, const tmp<block<Type,P>>&);

// template<Type1, Type2>, VS = VectorSpace, CS = CellSpace
//
//     block<arg1<Type1,Type2>> = block<Type1,P> arg2 block<Type2,P>
//     block<arg1<Type1,Type2>> = VS<Type1> arg2 block<Type2,P>
//     block<arg1<Type1,Type2>> = block<Type1,P> arg2 VS<Type2>
//     block<arg1<Type1,Type2>> = CS<Type1> arg2 block<Type2,P>
//     block<arg1<Type1,Type2>> = block<Type1,P> arg2 CS<Type2>
//
// Note: this does not define
//
//     block<arg1<Type1,Type2>> = Type1 arg2 block<Type2,P>
//     block<arg1<Type1,Type2>> = block<Type1,P> arg2 Type2
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
#include "labelBlock.H"
#include "labelVectorBlock.H"
#include "vectorBlock.H"
#include "tensorBlock.H"
#include "symmTensorBlock.H"
#include "faceScalarBlock.H"
#include "vertexScalarBlock.H"
#include "faceVectorBlock.H"
#include "vertexVectorBlock.H"
#include "diagStencilBlock.H"
#include "stencilBlock.H"
