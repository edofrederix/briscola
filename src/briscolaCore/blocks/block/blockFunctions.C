#include "PstreamReduceOps.H"
#include "blockReuseFunctions.H"

#define TEMPLATE template<class Type>
#include "blockFunctionsM.C"

// Scalar return type must be deduced because of cell space
#define SCALARPRODTYPE typename scalarProduct<Type,Type>::type

namespace Foam
{

namespace briscola
{

template<class Type>
void component
(
    block<typename block<Type>::cmptType>& res,
    const block<Type>& f,
    const direction d
)
{
    forAllBlockLinear(res, i)
        res(i) = f(i).component(d);
}

template<class Type>
void T(block<Type>& res, const block<Type>& f)
{
    forAllBlockLinear(res, i)
        res(i) = Foam::T(f(i));
}

template<class Type>
void sqr
(
    block<typename outerProduct<Type, Type>::type>& res,
    const block<Type>& vf
)
{
    forAllBlockLinear(res, i)
        res(i) = Foam::sqr(vf(i));
}

template<class Type>
tmp<block<typename outerProduct<Type, Type>::type>>
sqr(const block<Type>& f)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<block<outerProductType>> tRes
    (
        new block<outerProductType>(f.size())
    );
    sqr(tRes.ref(), f);
    return tRes;
}

template<class Type>
tmp<block<typename outerProduct<Type, Type>::type>>
sqr(const tmp<block<Type>>& tf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<block<outerProductType>> tRes =
        reuseTmp<outerProductType, Type>::New(tf);
    sqr(tRes.ref(), tf());
    tf.clear();
    return tRes;
}

template<class Type>
void mag(block<SCALARPRODTYPE>& res, const block<Type>& f)
{
    forAllBlockLinear(res, i)
        res(i) = Foam::mag(f(i));
}

template<class Type>
tmp<block<SCALARPRODTYPE>> mag(const block<Type>& f)
{
    tmp<block<SCALARPRODTYPE>> tRes(new block<SCALARPRODTYPE>(f.shape()));
    mag(tRes.ref(), f);
    return tRes;
}

template<class Type>
tmp<block<SCALARPRODTYPE>> mag(const tmp<block<Type>>& tf)
{
    tmp<block<SCALARPRODTYPE>> tRes = reuseTmp<SCALARPRODTYPE, Type>::New(tf);
    mag(tRes.ref(), tf());
    if (tf.isTmp())
        tf.clear();
    return tRes;
}

template<class Type>
void cmptMag(block<Type>& res, const block<Type>& f)
{
    forAllBlockLinear(res, i)
        res(i) = Foam::cmptMag(f(i));
}

template<class Type>
tmp<block<Type>> cmptMag(const block<Type>& f)
{
    tmp<block<Type>> tRes(new block<Type>(f.shape()));
    cmptMag(tRes.ref(), f);
    return tRes;
}

template<class Type>
tmp<block<Type>> cmptMag(const tmp<block<Type>>& tf)
{
    tmp<block<Type>> tRes = New(tf);
    cmptMag(tRes.ref(), tf());
    if (tf.isTmp())
        tf.clear();
    return tRes;
}

template<class Type>
void cmptSqr(block<Type>& res, const block<Type>& f)
{
    forAllBlockLinear(res, i)
        res(i) = Foam::cmptSqr(f(i));
}

template<class Type>
tmp<block<Type>> cmptSqr(const block<Type>& f)
{
    tmp<block<Type>> tRes(new block<Type>(f.shape()));
    cmptSqr(tRes.ref(), f);
    return tRes;
}

template<class Type>
tmp<block<Type>> cmptSqr(const tmp<block<Type>>& tf)
{
    tmp<block<Type>> tRes = New(tf);
    cmptSqr(tRes.ref(), tf());
    if (tf.isTmp())
        tf.clear();
    return tRes;
}

#define TMP_UNARY_FUNCTION(ReturnType, Func)                                   \
                                                                               \
template<class Type>                                                           \
ReturnType Func(const tmp<block<Type>>& tf1)                                   \
{                                                                              \
    ReturnType res = Func(tf1());                                              \
    if (tf1.isTmp()) tf1.clear();                                              \
    return res;                                                                \
}

template<class Type>
Type max(const block<Type>& f)
{
    if (f.size())
    {
        Type Max(f(0));

        forAllBlockLinear(f, i)
            if (f(i) > Max)
                Max = f(i);

        return Max;
    }
    else
    {
        return pTraits<Type>::min;
    }
}

TMP_UNARY_FUNCTION(Type, max)

template<class Type>
Type min(const block<Type>& f)
{
    if (f.size())
    {
        Type Min(f(0));

        forAllBlockLinear(f, i)
            if (f(i) < Min)
                Min = f(i);

        return Min;
    }
    else
    {
        return pTraits<Type>::max;
    }
}

TMP_UNARY_FUNCTION(Type, min)

template<class Type>
Type sum(const block<Type>& f)
{
    if (f.size())
    {
        Type Sum = Zero;

        forAllBlockLinear(f, i)
            Sum += f(i);

        return Sum;
    }
    else
    {
        return Zero;
    }
}

TMP_UNARY_FUNCTION(Type, sum)

template<class Type>
Type average(const block<Type>& f)
{
    if (f.size())
    {
        Type avrg = sum(f)/f.size();

        return avrg;
    }
    else
    {
        return Zero;
    }
}

TMP_UNARY_FUNCTION(Type, average)

#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                       \
                                                                               \
template<class Type>                                                           \
ReturnType gFunc                                                               \
(                                                                              \
    const block<Type>& f,                                                      \
    const label comm                                                           \
)                                                                              \
{                                                                              \
    ReturnType res = Func(f);                                                  \
    reduce(res, rFunc##Op<ReturnType>(), Pstream::msgType(), comm);            \
    return res;                                                                \
}                                                                              \
                                                                               \
template<class Type>                                                           \
ReturnType gFunc                                                               \
(                                                                              \
    const tmp<block<Type>>& tf,                                                \
    const label comm                                                           \
)                                                                              \
{                                                                              \
    return gFunc(tf());                                                        \
}

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)

#undef G_UNARY_FUNCTION

template<class Type>
Type gAverage
(
    const block<Type>& f,
    const label comm
)
{
    label n = f.size();
    Type s = sum(f);
    sumReduce(s, n, Pstream::msgType(), comm);

    if (n > 0)
    {
        return s/n;
    }
    else
    {
        return Zero;
    }
}

template<class Type>
Type gAverage
(
    const tmp<block<Type>>& tf,
    const label comm
)
{
    return gAverage(tf(),comm);
}

#undef TMP_UNARY_FUNCTION

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

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type1, Type2>::type>& res,                          \
    const block<Type1>& f1,                                                    \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    checkBlocks(res,f1,f2,#OpFunc);                                            \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = f1(i) Op f2(i);                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const block<Type1>& f1, const block<Type2>& f2)                    \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType>> tRes(new block<productType>(f1.shape()));          \
    OpFunc(tRes.ref(), f1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const block<Type1>& f1, const tmp<block<Type2>>& tf2)              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType>> tRes = reuseTmp<productType, Type2>::New(tf2);     \
    OpFunc(tRes.ref(), f1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<block<Type1>>& tf1, const block<Type2>& f2)              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType>> tRes = reuseTmp<productType, Type1>::New(tf1);     \
    OpFunc(tRes.ref(), tf1(), f2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<block<Type1>>& tf1, const tmp<block<Type2>>& tf2)        \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType>> tRes =                                             \
        reuseTmpTmp<productType, Type1, Type1, Type2>::New(tf1, tf2);          \
    OpFunc(tRes.ref(), tf1(), tf2());                                          \
    if (tf1.isTmp()) tf1.clear();                                              \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
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
    block<typename product<Type, Form>::type>& res,                            \
    const block<Type>& f1,                                                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    checkBlocks(res,f1,#OpFunc);                                               \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = f1(i) Op static_cast<const Form&>(vs);                        \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op(const block<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs)     \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType>> tRes(new block<productType>(f1.shape()));          \
    OpFunc(tRes.ref(), f1, static_cast<const Form&>(vs));                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const tmp<block<Type>>& tf1,                                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType>> tRes = reuseTmp<productType, Type>::New(tf1);      \
    OpFunc(tRes.ref(), tf1(), static_cast<const Form&>(vs));                   \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
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
    block<typename product<Form, Type>::type>& res,                            \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const block<Type>& f1                                                      \
)                                                                              \
{                                                                              \
    checkBlocks(res,f1,#OpFunc);                                               \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = static_cast<const Form&>(vs) Op f1(i);                        \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const block<Type>& f1)     \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType>> tRes(new block<productType>(f1.shape()));          \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), f1);                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<block<Type>>& tf1                                                \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType>> tRes = reuseTmp<productType, Type>::New(tf1);      \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), tf1());                   \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
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
    block<typename product<Type, Form>::type>& res,                            \
    const block<Type>& f1,                                                     \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    checkBlocks(res,f1,#OpFunc);                                               \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = f1(i) Op static_cast<const Form&>(vs);                        \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op(const block<Type>& f1, const CellSpace<Form,Cmpt,nCmpt>& vs)       \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType>> tRes(new block<productType>(f1.shape()));          \
    OpFunc(tRes.ref(), f1, static_cast<const Form&>(vs));                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const tmp<block<Type>>& tf1,                                               \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType>> tRes = reuseTmp<productType, Type>::New(tf1);      \
    OpFunc(tRes.ref(), tf1(), static_cast<const Form&>(vs));                   \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
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
    block<typename product<Form, Type>::type>& res,                            \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const block<Type>& f1                                                      \
)                                                                              \
{                                                                              \
    checkBlocks(res,f1,#OpFunc);                                               \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = static_cast<const Form&>(vs) Op f1(i);                        \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op(const CellSpace<Form,Cmpt,nCmpt>& vs, const block<Type>& f1)       \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType>> tRes(new block<productType>(f1.shape()));          \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), f1);                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const tmp<block<Type>>& tf1                                                \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType>> tRes = reuseTmp<productType, Type>::New(tf1);      \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), tf1());                   \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}

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
