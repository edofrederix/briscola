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
    typedef typename block<Type>::cmptType cmptType;
    BFOR_ALL_F_OP_F_FUNC_S
    (
        cmptType, res, =, Type, f, .component, const direction, d
    )
}


template<class Type>
void T(block<Type>& res, const block<Type>& f)
{
    BFOR_ALL_F_OP_F_FUNC(Type, res, =, Type, f, T)
}


template<class Type, direction r>
void pow
(
    block<typename powProduct<Type, r>::type>& res,
    const block<Type>& vf
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    BFOR_ALL_F_OP_FUNC_F_S
    (
        powProductType, res, =, pow, Type, vf, powProductType,
        pTraits<powProductType>::zero
    )
}

template<class Type, direction r>
tmp<block<typename powProduct<Type, r>::type>>
pow
(
    const block<Type>& f,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<block<powProductType>> tRes
    (
        new block<powProductType>(f.size())
    );
    pow<Type, r>(tRes.ref(), f);
    return tRes;
}

template<class Type, direction r>
tmp<block<typename powProduct<Type, r>::type>>
pow
(
    const tmp<block<Type>>& tf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<block<powProductType>> tRes = reuseTmp<powProductType, Type>::New(tf);
    pow<Type, r>(tRes.ref(), tf());
    tf.clear();
    return tRes;
}


template<class Type>
void sqr
(
    block<typename outerProduct<Type, Type>::type>& res,
    const block<Type>& vf
)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    BFOR_ALL_F_OP_FUNC_F(outerProductType, res, =, sqr, Type, vf)
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
void magSqr(block<scalar>& res, const block<Type>& f)
{
    BFOR_ALL_F_OP_FUNC_F(scalar, res, =, magSqr, Type, f)
}

template<class Type>
tmp<block<scalar>> magSqr(const block<Type>& f)
{
    tmp<block<scalar>> tRes(new block<scalar>(f.size()));
    magSqr(tRes.ref(), f);
    return tRes;
}

template<class Type>
tmp<block<scalar>> magSqr(const tmp<block<Type>>& tf)
{
    tmp<block<scalar>> tRes = reuseTmp<scalar, Type>::New(tf);
    magSqr(tRes.ref(), tf());
    tf.clear();
    return tRes;
}

template<class Type>
void mag(block<SCALARPRODTYPE>& res, const block<Type>& f)
{
    typedef SCALARPRODTYPE ScalarType;
    BFOR_ALL_F_OP_FUNC_F(ScalarType, res, =, ::Foam::mag, Type, f)
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
void cmptMax(block<typename block<Type>::cmptType>& res, const block<Type>& f)
{
    typedef typename block<Type>::cmptType cmptType;
    BFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptMax, Type, f)
}

template<class Type>
tmp<block<typename block<Type>::cmptType>> cmptMax(const block<Type>& f)
{
    typedef typename block<Type>::cmptType cmptType;
    tmp<block<cmptType>> tRes(new block<cmptType>(f.shape()));
    cmptMax(tRes.ref(), f);
    return tRes;
}

template<class Type>
tmp<block<typename block<Type>::cmptType>> cmptMax(const tmp<block<Type>>& tf)
{
    typedef typename block<Type>::cmptType cmptType;
    tmp<block<cmptType>> tRes = reuseTmp<cmptType, Type>::New(tf);
    cmptMax(tRes.ref(), tf());
    if (tf.isTmp())
        tf.clear();
    return tRes;
}


template<class Type>
void cmptMin(block<typename block<Type>::cmptType>& res, const block<Type>& f)
{
    typedef typename block<Type>::cmptType cmptType;
    BFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptMin, Type, f)
}

template<class Type>
tmp<block<typename block<Type>::cmptType>> cmptMin(const block<Type>& f)
{
    typedef typename block<Type>::cmptType cmptType;
    tmp<block<cmptType>> tRes(new block<cmptType>(f.shape()));
    cmptMin(tRes.ref(), f);
    return tRes;
}

template<class Type>
tmp<block<typename block<Type>::cmptType>> cmptMin(const tmp<block<Type>>& tf)
{
    typedef typename block<Type>::cmptType cmptType;
    tmp<block<cmptType>> tRes = reuseTmp<cmptType, Type>::New(tf);
    cmptMin(tRes.ref(), tf());
    if (tf.isTmp())
        tf.clear();
    return tRes;
}

template<class Type>
void cmptAv(block<typename block<Type>::cmptType>& res, const block<Type>& f)
{
    typedef typename block<Type>::cmptType cmptType;
    BFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptAv, Type, f)
}

template<class Type>
tmp<block<typename block<Type>::cmptType>> cmptAv(const block<Type>& f)
{
    typedef typename block<Type>::cmptType cmptType;
    tmp<block<cmptType>> tRes(new block<cmptType>(f.shape()));
    cmptAv(tRes.ref(), f);
    return tRes;
}

template<class Type>
tmp<block<typename block<Type>::cmptType>> cmptAv(const tmp<block<Type>>& tf)
{
    typedef typename block<Type>::cmptType cmptType;
    tmp<block<cmptType>> tRes = reuseTmp<cmptType, Type>::New(tf);
    cmptAv(tRes.ref(), tf());
    if (tf.isTmp())
        tf.clear();
    return tRes;
}

template<class Type>
void cmptMag(block<Type>& res, const block<Type>& f)
{
    BFOR_ALL_F_OP_FUNC_F(Type, res, =, cmptMag, Type, f)
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
    BFOR_ALL_F_OP_FUNC_F(Type, res, =, ::Foam::cmptSqr, Type, f)
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

template<class Type>
void cmptSqrt(block<Type>& res, const block<Type>& f)
{
    BFOR_ALL_F_OP_FUNC_F(Type, res, =, ::Foam::cmptSqrt, Type, f)
}

template<class Type>
tmp<block<Type>> cmptSqrt(const block<Type>& f)
{
    tmp<block<Type>> tRes(new block<Type>(f.shape()));
    cmptSqrt(tRes.ref(), f);
    return tRes;
}

template<class Type>
tmp<block<Type>> cmptSqrt(const tmp<block<Type>>& tf)
{
    tmp<block<Type>> tRes = New(tf);
    cmptSqrt(tRes.ref(), tf());
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
        BFOR_ALL_S_OP_FUNC_F_S(Type, Max, =, ::Foam::max, Type, f, Type, Max)
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
        BFOR_ALL_S_OP_FUNC_F_S(Type, Min, =, ::Foam::min, Type, f, Type, Min)
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
        BFOR_ALL_S_OP_F(Type, Sum, +=, Type, f)
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

template<class Type>
Type maxMagSqr(const block<Type>& f)
{
    if (f.size())
    {
        Type Max(f(0));
        BFOR_ALL_S_OP_FUNC_F_S
        (
            Type,
            Max,
            =,
            ::Foam::maxMagSqrOp<Type>(),
            Type,
            f,
            Type,
            Max
        )
        return Max;
    }
    else
    {
        return Zero;
    }
}

TMP_UNARY_FUNCTION(Type, maxMagSqr)

template<class Type>
Type minMagSqr(const block<Type>& f)
{
    if (f.size())
    {
        Type Min(f(0));
        BFOR_ALL_S_OP_FUNC_F_S
        (
            Type,
            Min,
            =,
            ::Foam::minMagSqrOp<Type>(),
            Type,
            f,
            Type,
            Min
        )
        return Min;
    }
    else
    {
        return pTraits<Type>::rootMax;
    }
}

TMP_UNARY_FUNCTION(Type, minMagSqr)

template<class Type>
SCALARPRODTYPE sumProd(const block<Type>& f1, const block<Type>& f2)
{
    if (f1.size() && (f1.size() == f2.size()))
    {
        typedef SCALARPRODTYPE ScalarType;
        ScalarType SumProd = Zero;
        BFOR_ALL_S_OP_F_OP_F(ScalarType, SumProd, +=, Type, f1, &&, Type, f2)
        return SumProd;
    }
    else
    {
        return Zero;
    }
}

template<class Type>
SCALARPRODTYPE
sumProd(const tmp<block<Type>>& tf1, const block<Type>& f2)
{
    return sumProd(tf1(),f2);
}

template<class Type>
SCALARPRODTYPE
sumProd(const block<Type>& f1, const tmp<block<Type>>& tf2)
{
    return sumProd(f1,tf2());
}

template<class Type>
SCALARPRODTYPE
sumProd(const tmp<block<Type>>& tf1, const tmp<block<Type>>& tf2)
{
    return sumProd(tf1(),tf2());
}

template<class Type>
Type sumCmptProd(const block<Type>& f1, const block<Type>& f2)
{
    if (f1.size() && (f1.size() == f2.size()))
    {
        Type SumProd = Zero;

        BFOR_ALL_S_OP_FUNC_F_F
        (
            Type,
            SumProd,
            +=,
            cmptMultiply,
            Type,
            f1,
            Type,
            f2
        )

        return SumProd;
    }
    else
    {
        return Zero;
    }
}

template<class Type>
Type sumCmptProd(const tmp<block<Type>>& tf1, const block<Type>& f2)
{
    return sumCmptProd(tf1(),f2);
}

template<class Type>
Type sumCmptProd(const block<Type>& f1, const tmp<block<Type>>& tf2)
{
    return sumCmptProd(f1,tf2());
}

template<class Type>
Type sumCmptProd(const tmp<block<Type>>& tf1, const tmp<block<Type>>& tf2)
{
    return sumCmptProd(tf1(),tf2());
}

template<class Type>
typename powProduct<Type,2>::type sumSqr(const block<Type>& f)
{
    if (f.size())
    {
        typedef typename powProduct<Type,2>::type ProductType;
        ProductType SumSqr = Zero;
        BFOR_ALL_S_OP_FUNC_F(ProductType, SumSqr, +=, ::Foam::sqr, Type, f)
        return SumSqr;
    }
    else
    {
        return Zero;
    }
}

template<class Type>
typename powProduct<Type,2>::type sumSqr(const tmp<block<Type>>& tf)
{
    return sumSqr(tf());
}

template<class Type>
SCALARPRODTYPE sumMag(const block<Type>& f)
{
    if (f.size())
    {
        typedef SCALARPRODTYPE ScalarType;
        ScalarType SumMag = Zero;
        BFOR_ALL_S_OP_FUNC_F(scalar, SumMag, +=, ::Foam::mag, Type, f)
        return SumMag;
    }
    else
    {
        return Zero;
    }
}

template<class Type>
SCALARPRODTYPE sumMag(const tmp<block<Type>>& tf1)
{
    SCALARPRODTYPE res = sumMag(tf1());
    if (tf1.isTmp())
        tf1.clear();
    return res;
}

template<class Type>
Type sumCmptMag(const block<Type>& f)
{
    if (f.size())
    {
        Type SumMag = Zero;
        BFOR_ALL_S_OP_FUNC_F(scalar, SumMag, +=, cmptMag, Type, f)
        return SumMag;
    }
    else
    {
        return Zero;
    }
}

TMP_UNARY_FUNCTION(Type, sumCmptMag)


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
G_UNARY_FUNCTION(Type, gMaxMagSqr, maxMagSqr, maxMagSqr)
G_UNARY_FUNCTION(Type, gMinMagSqr, minMagSqr, minMagSqr)
G_UNARY_FUNCTION(SCALARPRODTYPE, gSumMag, sumMag, sum)
G_UNARY_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)

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

template<class Type>
typename powProduct<Type,2>::type gSumSqr
(
    const block<Type>& f,
    const label comm
)
{
    typename powProduct<Type,2>::type res = sumSqr(f);
    reduce
    (
        res,
        sumOp<typename powProduct<Type,2>::type>(),
        Pstream::msgType(),
        comm
    );
    return res;
}

template<class Type>
typename powProduct<Type,2>::type gSumSqr
(
    const tmp<block<Type>>& tf,
    const label comm
)
{
    return gSumSqr(tf());
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
