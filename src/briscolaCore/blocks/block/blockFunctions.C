#include "PstreamReduceOps.H"
#include "blockReuseFunctions.H"

#define TEMPLATE template<class Type, int P>
#include "blockFunctionsM.C"

// Scalar return type must be deduced because of cell space
#define SCALARPRODTYPE typename scalarProduct<Type,Type>::type

namespace Foam
{

namespace briscola
{

template<class Type, int P>
void component
(
    block<typename block<Type,P>::cmptType, P>& res,
    const block<Type,P>& f,
    const direction d
)
{
    forAllBlockLinear(res, i)
        res(i) = f(i).component(d);
}

template<class Type, int P>
void T(block<Type,P>& res, const block<Type,P>& f)
{
    forAllBlockLinear(res, i)
        res(i) = Foam::T(f(i));
}

template<class Type, int P>
void sqr
(
    block<typename outerProduct<Type, Type>::type, P>& res,
    const block<Type,P>& vf
)
{
    forAllBlockLinear(res, i)
        res(i) = Foam::sqr(vf(i));
}

template<class Type, int P>
tmp<block<typename outerProduct<Type, Type>::type, P>>
sqr(const block<Type,P>& f)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<block<outerProductType, P>> tRes
    (
        new block<outerProductType, P>(f.size())
    );
    sqr(tRes.ref(), f);
    return tRes;
}

template<class Type, int P>
tmp<block<typename outerProduct<Type, Type>::type, P>>
sqr(const tmp<block<Type,P>>& tf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<block<outerProductType, P>> tRes =
        reuseTmp<outerProductType, Type, P>::New(tf);
    sqr(tRes.ref(), tf());
    tf.clear();
    return tRes;
}

template<class Type, int P>
void mag(block<SCALARPRODTYPE, P>& res, const block<Type,P>& f)
{
    forAllBlockLinear(res, i)
        res(i) = Foam::mag(f(i));
}

template<class Type, int P>
tmp<block<SCALARPRODTYPE, P>> mag(const block<Type,P>& f)
{
    tmp<block<SCALARPRODTYPE, P>> tRes(new block<SCALARPRODTYPE, P>(f.shape()));
    mag(tRes.ref(), f);
    return tRes;
}

template<class Type, int P>
tmp<block<SCALARPRODTYPE, P>> mag(const tmp<block<Type,P>>& tf)
{
    tmp<block<SCALARPRODTYPE, P>> tRes =
        reuseTmp<SCALARPRODTYPE, Type, P>::New(tf);
    mag(tRes.ref(), tf());
    if (tf.isTmp())
        tf.clear();
    return tRes;
}

template<class Type, int P>
void cmptMag(block<Type,P>& res, const block<Type,P>& f)
{
    forAllBlockLinear(res, i)
        res(i) = Foam::cmptMag(f(i));
}

template<class Type, int P>
tmp<block<Type,P>> cmptMag(const block<Type,P>& f)
{
    tmp<block<Type,P>> tRes(new block<Type,P>(f.shape()));
    cmptMag(tRes.ref(), f);
    return tRes;
}

template<class Type, int P>
tmp<block<Type,P>> cmptMag(const tmp<block<Type,P>>& tf)
{
    tmp<block<Type,P>> tRes = New(tf);
    cmptMag(tRes.ref(), tf());
    if (tf.isTmp())
        tf.clear();
    return tRes;
}

template<class Type, int P>
void cmptSqr(block<Type,P>& res, const block<Type,P>& f)
{
    forAllBlockLinear(res, i)
        res(i) = Foam::cmptSqr(f(i));
}

template<class Type, int P>
tmp<block<Type,P>> cmptSqr(const block<Type,P>& f)
{
    tmp<block<Type,P>> tRes(new block<Type,P>(f.shape()));
    cmptSqr(tRes.ref(), f);
    return tRes;
}

template<class Type, int P>
tmp<block<Type,P>> cmptSqr(const tmp<block<Type,P>>& tf)
{
    tmp<block<Type,P>> tRes = New(tf);
    cmptSqr(tRes.ref(), tf());
    if (tf.isTmp())
        tf.clear();
    return tRes;
}

#define TMP_UNARY_FUNCTION(ReturnType, Func)                                   \
                                                                               \
template<class Type, int P>                                                    \
ReturnType Func(const tmp<block<Type,P>>& tf1)                                 \
{                                                                              \
    ReturnType res = Func(tf1());                                              \
    if (tf1.isTmp()) tf1.clear();                                              \
    return res;                                                                \
}

template<class Type, int P>
Type max(const block<Type,P>& f)
{
    if (f.size())
    {
        Type Max(f(0,0,0));

        forAllBlock(f, i, j, k)
            if (f(i,j,k) > Max)
                Max = f(i,j,k);

        return Max;
    }
    else
    {
        return pTraits<Type>::min;
    }
}

TMP_UNARY_FUNCTION(Type, max)

template<class Type, int P>
Type min(const block<Type,P>& f)
{
    if (f.size())
    {
        Type Min(f(0,0,0));

        forAllBlock(f, i, j, k)
            if (f(i,j,k) < Min)
                Min = f(i,j,k);

        return Min;
    }
    else
    {
        return pTraits<Type>::max;
    }
}

TMP_UNARY_FUNCTION(Type, min)

template<class Type, int P>
Type sum(const block<Type,P>& f)
{
    if (f.size())
    {
        Type Sum = Zero;

        forAllBlock(f, i, j, k)
            Sum += f(i,j,k);

        return Sum;
    }
    else
    {
        return Zero;
    }
}

TMP_UNARY_FUNCTION(Type, sum)

template<class Type, int P>
Type average(const block<Type,P>& f)
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
template<class Type, int P>                                                    \
ReturnType gFunc                                                               \
(                                                                              \
    const block<Type,P>& f,                                                    \
    const label comm                                                           \
)                                                                              \
{                                                                              \
    ReturnType res = Func(f);                                                  \
    reduce(res, rFunc##Op<ReturnType>(), Pstream::msgType(), comm);            \
    return res;                                                                \
}                                                                              \
                                                                               \
template<class Type, int P>                                                    \
ReturnType gFunc                                                               \
(                                                                              \
    const tmp<block<Type,P>>& tf,                                              \
    const label comm                                                           \
)                                                                              \
{                                                                              \
    return gFunc(tf());                                                        \
}

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)

#undef G_UNARY_FUNCTION

template<class Type, int P>
Type gAverage
(
    const block<Type,P>& f,
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

template<class Type, int P>
Type gAverage
(
    const tmp<block<Type,P>>& tf,
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

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2, int P>                                      \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type1, Type2>::type, P>& res,                       \
    const block<Type1,P>& f1,                                                  \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type ReturnType;                   \
    checkBlocks<ReturnType,Type1,Type2,P>(res,f1,f2,#OpFunc);                  \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = f1(i) Op f2(i);                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, int P>                                      \
tmp<block<typename product<Type1, Type2>::type, P>>                            \
operator Op(const block<Type1,P>& f1, const block<Type2,P>& f2)                \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType, P>> tRes(new block<productType, P>(f1.shape()));    \
    OpFunc(tRes.ref(), f1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, int P>                                      \
tmp<block<typename product<Type1, Type2>::type, P>>                            \
operator Op(const block<Type1,P>& f1, const tmp<block<Type2,P>>& tf2)          \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType, P>> tRes =                                          \
        reuseTmp<productType, Type2, P>::New(tf2);                             \
    OpFunc(tRes.ref(), f1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, int P>                                      \
tmp<block<typename product<Type1, Type2>::type, P>>                            \
operator Op(const tmp<block<Type1,P>>& tf1, const block<Type2,P>& f2)          \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType,P>> tRes = reuseTmp<productType, Type1, P>::New(tf1);\
    OpFunc(tRes.ref(), tf1(), f2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, int P>                                      \
tmp<block<typename product<Type1, Type2>::type, P>>                            \
operator Op(const tmp<block<Type1,P>>& tf1, const tmp<block<Type2,P>>& tf2)    \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType, P>> tRes =                                          \
        reuseTmpTmp<productType, Type1, Type1, Type2, P>::New(tf1, tf2);       \
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
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type, Form>::type, P>& res,                         \
    const block<Type,P>& f1,                                                   \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type ReturnType;                     \
    checkBlocks<ReturnType,Type,P>(res,f1,#OpFunc);                            \
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
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Type, Form>::type, P>>                              \
operator Op(const block<Type,P>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs)   \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType, P>> tRes(new block<productType, P>(f1.shape()));    \
    OpFunc(tRes.ref(), f1, static_cast<const Form&>(vs));                      \
    return tRes;                                                               \
}                                                                              \
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
operator Op                                                                    \
(                                                                              \
    const tmp<block<Type,P>>& tf1,                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType, P>> tRes =                                          \
        reuseTmp<productType, Type, P>::New(tf1);                              \
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
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Form, Type>::type, P>& res,                         \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const block<Type,P>& f1                                                    \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type ReturnType;                     \
    checkBlocks<ReturnType,Type,P>(res,f1,#OpFunc);                            \
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
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Form, Type>::type, P>>                              \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const block<Type,P>& f1)   \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType, P>> tRes(new block<productType, P>(f1.shape()));    \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), f1);                      \
    return tRes;                                                               \
}                                                                              \
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
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<block<Type,P>>& tf1                                              \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType, P>> tRes = reuseTmp<productType, Type, P>::New(tf1);\
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
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type, Form>::type, P>& res,                         \
    const block<Type,P>& f1,                                                   \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type ReturnType;                     \
    checkBlocks<ReturnType,Type,P>(res,f1,#OpFunc);                            \
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
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Type, Form>::type, P>>                              \
operator Op(const block<Type,P>& f1, const CellSpace<Form,Cmpt,nCmpt>& vs)     \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType,P>> tRes(new block<productType,P>(f1.shape()));      \
    OpFunc(tRes.ref(), f1, static_cast<const Form&>(vs));                      \
    return tRes;                                                               \
}                                                                              \
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
operator Op                                                                    \
(                                                                              \
    const tmp<block<Type,P>>& tf1,                                             \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType,P>> tRes = reuseTmp<productType, Type, P>::New(tf1); \
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
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Form, Type>::type, P>& res,                         \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const block<Type,P>& f1                                                    \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type ReturnType;                     \
    checkBlocks<ReturnType,Type,P>(res,f1,#OpFunc);                            \
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
    direction nCmpt,                                                           \
    int P                                                                      \
>                                                                              \
tmp<block<typename product<Form, Type>::type, P>>                              \
operator Op(const CellSpace<Form,Cmpt,nCmpt>& vs, const block<Type,P>& f1)     \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType,P>> tRes(new block<productType,P>(f1.shape()));      \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), f1);                      \
    return tRes;                                                               \
}                                                                              \
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
operator Op                                                                    \
(                                                                              \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const tmp<block<Type,P>>& tf1                                              \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType,P>> tRes = reuseTmp<productType, Type, P>::New(tf1); \
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
