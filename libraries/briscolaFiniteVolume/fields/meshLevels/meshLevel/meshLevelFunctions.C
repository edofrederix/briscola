#include "PstreamReduceOps.H"
#include "meshLevelReuseFunctions.H"

#define TEMPLATE template<class Type, class MeshType>
#include "meshLevelFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void mag(meshLevel<scalar,MeshType>& res, const meshLevel<Type,MeshType>& f)
{
    forAll(res, d)
        mag(res[d], f[d]);
}

template<class Type, class MeshType>
tmp<meshLevel<scalar,MeshType>> mag(const meshLevel<Type,MeshType>& f)
{
    tmp<meshLevel<scalar,MeshType>>
        tRes(new meshLevel<scalar,MeshType>(f.fvMsh(),f.levelNum()));
    mag(tRes.ref(), f);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshLevel<scalar,MeshType>> mag(const tmp<meshLevel<Type,MeshType>>& tf)
{
    tmp<meshLevel<scalar,MeshType>> tRes =
        reuseLevTmp<scalar,Type,MeshType>::New(tf);

    mag(tRes.ref(), tf());
    tf.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptMax
(
    meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>& res,
    const meshLevel<Type,MeshType>& f
)
{
    forAll(res, d)
        cmptMax(res[d], f[d]);
}

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptMax(const meshLevel<Type,MeshType>& f)
{
    typedef typename meshLevel<Type,MeshType>::cmptType cmptType;
    tmp<meshLevel<cmptType,MeshType>>
        tRes(new meshLevel<cmptType,MeshType>(f.fvMsh(),f.levelNum()));
    cmptMax(tRes.ref(), f);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptMax(const tmp<meshLevel<Type,MeshType>>& tf)
{
    typedef typename meshLevel<Type,MeshType>::cmptType cmptType;
    tmp<meshLevel<cmptType,MeshType>> tRes =
        reuseLevTmp<cmptType,Type,MeshType>::New(tf);

    cmptMax(tRes.ref(), tf());
    tf.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptMin
(
    meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>& res,
    const meshLevel<Type,MeshType>& f
)
{
    forAll(res, d)
        cmptMin(res[d], f[d]);
}

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptMin(const meshLevel<Type,MeshType>& f)
{
    typedef typename meshLevel<Type,MeshType>::cmptType cmptType;
    tmp<meshLevel<cmptType,MeshType>>
        tRes(new meshLevel<cmptType,MeshType>(f.fvMsh(),f.levelNum()));
    cmptMin(tRes.ref(), f);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptMin(const tmp<meshLevel<Type,MeshType>>& tf)
{
    typedef typename meshLevel<Type,MeshType>::cmptType cmptType;
    tmp<meshLevel<cmptType,MeshType>> tRes =
        reuseLevTmp<cmptType,Type,MeshType>::New(tf);

    cmptMin(tRes.ref(), tf());
    tf.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptAv
(
    meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>& res,
    const meshLevel<Type,MeshType>& f
)
{
    forAll(res, d)
        cmptAv(res[d],f[d]);
}

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptAv(const meshLevel<Type,MeshType>& f)
{
    typedef typename meshLevel<Type,MeshType>::cmptType cmptType;
    tmp<meshLevel<cmptType,MeshType>>
        tRes(new meshLevel<cmptType,MeshType>(f.fvMsh(),f.levelNum()));
    cmptAv(tRes.ref(), f);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
cmptAv(const tmp<meshLevel<Type,MeshType>>& tf)
{
    typedef typename meshLevel<Type,MeshType>::cmptType cmptType;
    tmp<meshLevel<cmptType,MeshType>> tRes =
        reuseLevTmp<cmptType,Type,MeshType>::New(tf);

    cmptAv(tRes.ref(), tf());
    tf.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptMag(meshLevel<Type,MeshType>& res, const meshLevel<Type,MeshType>& f)
{
    forAll(res, d)
        cmptMag(res[d],f[d]);
}

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptMag(const meshLevel<Type,MeshType>& f)
{
    tmp<meshLevel<Type,MeshType>>
        tRes(new meshLevel<Type,MeshType>(f.fvMsh(),f.levelNum()));
    cmptMag(tRes.ref(), f);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptMag(const tmp<meshLevel<Type,MeshType>>& tf)
{
    tmp<meshLevel<Type,MeshType>> tRes = New(tf);

    cmptMag(tRes.ref(), tf());
    tf.clear();
    return tRes;
}

template<class Type, class MeshType>
List<Type> max(const meshLevel<Type,MeshType>& f)
{
    List<Type> Max(f.size());

    forAll(Max, d)
    {
        Max[d] = max(f[d]);
    }

    return Max;
}

template<class Type, class MeshType>
List<Type> max(const tmp<meshLevel<Type,MeshType>>& tf)
{
    List<Type> ret(max(tf()));
    tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type> min(const meshLevel<Type,MeshType>& f)
{
    List<Type> Min(f.size());

    forAll(Min, d)
    {
        Min[d] = min(f[d]);
    }

    return Min;
}

template<class Type, class MeshType>
List<Type> min(const tmp<meshLevel<Type,MeshType>>& tf)
{
    List<Type> ret(min(tf()));
    tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type> sum(const meshLevel<Type,MeshType>& f)
{
    List<Type> Sum(f.size());

    forAll(Sum, d)
    {
        Sum[d] = sum(f[d]);
    }

    return Sum;
}

template<class Type, class MeshType>
List<Type> sum(const tmp<meshLevel<Type,MeshType>>& tf)
{
    List<Type> ret(sum(tf()));
    tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type> average(const meshLevel<Type,MeshType>& f)
{
    // Only for interior cells

    List<Type> Sum(sum(f));

    forAll(f, d)
        Sum[d] /= f[d].size();

    return Sum;
}

template<class Type, class MeshType>
List<Type> average(const tmp<meshLevel<Type,MeshType>>& tf)
{
    List<Type> ret(average(tf()));
    tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<scalar> maxMagSqr(const meshLevel<Type,MeshType>& f)
{
    List<scalar> Max(f.size());

    forAll(Max, d)
    {
        Max[d] = maxMagSqr(f[d]);
    }

    return Max;
}

template<class Type, class MeshType>
List<scalar> maxMagSqr(const tmp<meshLevel<Type,MeshType>>& tf)
{
    List<scalar> ret(maxMagSqr(tf()));
    tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<scalar> minMagSqr(const meshLevel<Type,MeshType>& f)
{
    List<scalar> Min(f.size());

    forAll(Min, d)
    {
        Min[d] = minMagSqr(f[d]);
    }

    return Min;
}

template<class Type, class MeshType>
List<scalar> minMagSqr(const tmp<meshLevel<Type,MeshType>>& tf)
{
    List<scalar> ret(minMagSqr(tf()));
    tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<scalar> sumProd
(
    const meshLevel<Type,MeshType>& f1,
    const meshLevel<Type,MeshType>& f2
)
{
    List<scalar> SumProd(f1.size());

    forAll(SumProd, d)
    {
        SumProd[d] = sumProd(f1[d], f2[d]);
    }

    return SumProd;
}

template<class Type, class MeshType>
List<scalar> sumProd
(
    const tmp<meshLevel<Type,MeshType>>& tf1,
    const meshLevel<Type,MeshType>& f2
)
{
    List<scalar> ret(sumProd(tf1(),f2));
    tf1.clear();
    return ret;
}

template<class Type, class MeshType>
List<scalar> sumProd
(
    const meshLevel<Type,MeshType>& f1,
    const tmp<meshLevel<Type,MeshType>>& tf2
)
{
    List<scalar> ret(sumProd(f1,tf2()));
    tf2.clear();
    return ret;
}

template<class Type, class MeshType>
List<scalar> sumProd
(
    const tmp<meshLevel<Type,MeshType>>& tf1,
    const tmp<meshLevel<Type,MeshType>>& tf2
)
{
    List<scalar> ret(sumProd(tf1(),tf2()));
    tf1.clear();
    tf2.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const meshLevel<Type,MeshType>& f1,
    const meshLevel<Type,MeshType>& f2
)
{
    List<Type> Sum(f1.size());

    forAll(Sum, d)
    {
        Sum[d] = sumCmptProd(f1[d],f2[d]);
    }

    return Sum;
}

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const tmp<meshLevel<Type,MeshType>>& tf1,
    const meshLevel<Type,MeshType>& f2
)
{
    List<Type> ret(sumCmptProd(tf1(),f2));
    tf1.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const meshLevel<Type,MeshType>& f1,
    const tmp<meshLevel<Type,MeshType>>& tf2
)
{
    List<Type> ret(sumCmptProd(f1,tf2()));
    tf2.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const tmp<meshLevel<Type,MeshType>>& tf1,
    const tmp<meshLevel<Type,MeshType>>& tf2
)
{
    List<Type> ret(sumCmptProd(tf1(),tf2()));
    tf1.clear();
    tf2.clear();
    return ret;
}

template<class Type, class MeshType>
List<scalar> sumMag(const meshLevel<Type,MeshType>& f)
{
    List<scalar> SumMag(f.size());

    forAll(SumMag, d)
    {
        SumMag[d] = sumMag(f[d]);
    }

    return SumMag;
}

template<class Type, class MeshType>
List<scalar> sumMag(const tmp<meshLevel<Type,MeshType>>& tf)
{
    List<scalar> ret(sumMag(tf()));
    tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type> sumCmptMag(const meshLevel<Type,MeshType>& f)
{
    List<Type> Sum(f.size());

    forAll(Sum, d)
    {
        Sum[d] = sumCmptMag(f[d]);
    }

    return Sum;
}

template<class Type, class MeshType>
List<Type> sumCmptMag(const tmp<meshLevel<Type,MeshType>>& tf)
{
    List<Type> ret(sumCmptMag(tf()));
    tf.clear();
    return ret;
}


#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                        \
                                                                                \
template<class Type, class MeshType>                                            \
List<ReturnType> gFunc                                                          \
(                                                                               \
    const meshLevel<Type,MeshType>& f,                                          \
    const label comm                                                            \
)                                                                               \
{                                                                               \
    List<ReturnType> res(Func(f));                                              \
    forAll(res, d)                                                              \
        reduce(res[d], rFunc##Op<ReturnType>(), Pstream::msgType(), comm);      \
    return res;                                                                 \
}                                                                               \
                                                                                \
template<class Type, class MeshType>                                            \
List<ReturnType> gFunc                                                          \
(                                                                               \
    const tmp<meshLevel<Type,MeshType>>& tf,                                    \
    const label comm                                                            \
)                                                                               \
{                                                                               \
    List<ReturnType> ret(gFunc(tf()));                                          \
    tf.clear();                                                                 \
    return ret;                                                                 \
}

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(scalar, gMaxMagSqr, maxMagSqr, maxMagSqr)
G_UNARY_FUNCTION(scalar, gMinMagSqr, minMagSqr, minMagSqr)
G_UNARY_FUNCTION(scalar, gSumMag, sumMag, sum)
G_UNARY_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)

#undef G_UNARY_FUNCTION

template<class Type, class MeshType>
List<Type> gAverage
(
    const meshLevel<Type,MeshType>& f,
    const label comm
)
{
    List<Type> Sum(sum(f));

    List<Type> Avrg
    (
        f.size(),
        Zero
    );

    forAll(f, d)
    {
        label n(cmptProduct(f[d].N()));

        sumReduce(Sum[d], n, Pstream::msgType(), comm);

        if (n > 0)
        {
            Avrg[d] = Sum[d]/n;
        }
    }

    return Avrg;
}

template<class Type, class MeshType>
List<Type> gAverage
(
    const tmp<meshLevel<Type,MeshType>>& tf,
    const label comm
)
{
    List<Type> ret(gAverage(tf()));
    tf.clear();
    return ret;
}

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
