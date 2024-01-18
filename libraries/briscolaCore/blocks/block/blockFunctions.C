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
void mag(block<SCALARPRODTYPE>& res, const block<Type>& f)
{
    forAllBlockLinear(f, i)
    {
        res(i) = Foam::mag(f(i));
    }
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
    tf.clear();
    return tRes;
}

template<class Type>
void cmptMax(block<typename block<Type>::cmptType>& res, const block<Type>& f)
{
    forAllBlockLinear(res, i)
    {
        res(i) = Foam::cmptMax(f(i));
    }
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
    tf.clear();
    return tRes;
}


template<class Type>
void cmptMin(block<typename block<Type>::cmptType>& res, const block<Type>& f)
{
    forAllBlockLinear(res, i)
    {
        res(i) = Foam::cmptMin(f(i));
    }
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
    tf.clear();
    return tRes;
}

template<class Type>
void cmptAv(block<typename block<Type>::cmptType>& res, const block<Type>& f)
{
    forAllBlockLinear(res, i)
    {
        res(i) = Foam::cmptAv(f(i));
    }
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
    tf.clear();
    return tRes;
}

template<class Type>
void cmptMag(block<Type>& res, const block<Type>& f)
{
    forAllBlockLinear(res, i)
    {
        res(i) = Foam::cmptMag(f(i));
    }
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
    tf.clear();
    return tRes;
}

#define TMP_UNARY_FUNCTION(ReturnType, Func)                                   \
                                                                               \
template<class Type>                                                           \
ReturnType Func(const tmp<block<Type>>& tf1)                                   \
{                                                                              \
    ReturnType res = Func(tf1());                                              \
    tf1.clear();                                                               \
    return res;                                                                \
}

template<class Type>
Type max(const block<Type>& f)
{
    if (f.size())
    {
        Type Max(f(0));

        forAllBlockLinear(f, i)
        {
            if (f(i) > Max)
            {
                Max = f(i);
            }
        }

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
        {
            if (f(i) < Min)
            {
                Min = f(i);
            }
        }

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
        {
            Sum += f(i);
        }

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
SCALARPRODTYPE maxMagSqr(const block<Type>& f)
{
    if (f.size())
    {
        SCALARPRODTYPE Max(Foam::magSqr(f(0)));

        forAllBlockLinear(f, i)
        {
            if (Foam::magSqr(f(i)) > Max)
            {
                Max = Foam::magSqr(f(i));
            }
        }

        return Max;
    }
    else
    {
        return Zero;
    }
}

template<class Type>
SCALARPRODTYPE maxMagSqr(const tmp<block<Type>>& tf1)
{
    SCALARPRODTYPE res = maxMagSqr(tf1());
    tf1.clear();
    return res;
}

template<class Type>
SCALARPRODTYPE minMagSqr(const block<Type>& f)
{
    if (f.size())
    {
        SCALARPRODTYPE Min(Foam::magSqr(f(0)));

        forAllBlockLinear(f, i)
        {
            if (Foam::magSqr(f(i)) < Min)
            {
                Min = Foam::magSqr(f(i));
            }
        }

        return Min;
    }
    else
    {
        return pTraits<SCALARPRODTYPE>::rootMax;
    }
}

template<class Type>
SCALARPRODTYPE minMagSqr(const tmp<block<Type>>& tf1)
{
    SCALARPRODTYPE res = minMagSqr(tf1());
    tf1.clear();
    return res;
}

template<class Type>
SCALARPRODTYPE sumProd(const block<Type>& f1, const block<Type>& f2)
{
    if (f1.size() && (f1.size() == f2.size()))
    {
        SCALARPRODTYPE SumProd = Zero;

        forAllBlockLinear(f1, i)
        {
            SumProd += (f1(i) && f2(i));
        }

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

        forAllBlockLinear(f1, i)
        {
            SumProd += Foam::cmptMultiply(f1(i),f2(i));
        }

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
        typename powProduct<Type,2>::type SumSqr = Zero;

        forAllBlockLinear(f, i)
        {
            SumSqr += Foam::sqr(f(i));
        }

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
        SCALARPRODTYPE SumMag = Zero;

        forAllBlockLinear(f, i)
        {
            SumMag += Foam::mag(f(i));
        }

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
    tf1.clear();
    return res;
}

template<class Type>
Type sumCmptMag(const block<Type>& f)
{
    if (f.size())
    {
        Type SumMag = Zero;

        forAllBlockLinear(f, i)
        {
            SumMag += Foam::cmptMag(f(i));
        }

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
G_UNARY_FUNCTION(SCALARPRODTYPE, gMaxMagSqr, maxMagSqr, maxMagSqr)
G_UNARY_FUNCTION(SCALARPRODTYPE, gMinMagSqr, minMagSqr, minMagSqr)
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
