#include "PstreamReduceOps.H"
#include "meshDirectionReuseFunctions.H"

#define TEMPLATE template<class Type, class MeshType>
#include "meshDirectionFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void mag
(
    meshDirection<scalar,MeshType>& res,
    const meshDirection<Type,MeshType>& D
)
{
    mag(res.B(), D.B());
}

template<class Type, class MeshType>
tmp<meshDirection<scalar,MeshType>>
mag(const meshDirection<Type,MeshType>& D)
{
    tmp<meshDirection<scalar,MeshType>> tRes
    (
        new meshDirection<scalar,MeshType>
        (
            D.fvMsh(),
            D.levelNum(),
            D.directionNum()
        )
    );

    mag(tRes.ref(), D);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshDirection<scalar,MeshType>>
mag(const tmp<meshDirection<Type,MeshType>>& tD)
{
    tmp<meshDirection<scalar,MeshType>> tRes =
        reuseDirTmp<scalar,Type,MeshType>::New(tD);

    mag(tRes.ref(), tD());
    tD.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptMax
(
    meshDirection
    <
        typename meshDirection<Type,MeshType>::cmptType,
        MeshType
    >& res,
    const meshDirection<Type,MeshType>& D
)
{
    cmptMax(res.B(), D.B());
}

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptMax(const meshDirection<Type,MeshType>& D)
{
    typedef typename meshDirection<Type,MeshType>::cmptType cmptType;

    tmp<meshDirection<cmptType,MeshType>> tRes
    (
        new meshDirection<cmptType,MeshType>
        (
            D.fvMsh(),
            D.levelNum(),
            D.directionNum()
        )
    );

    cmptMax(tRes.ref(), D);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptMax(const tmp<meshDirection<Type,MeshType>>& tD)
{
    typedef typename meshDirection<Type,MeshType>::cmptType cmptType;
    tmp<meshDirection<cmptType,MeshType>> tRes =
        reuseDirTmp<cmptType,Type,MeshType>::New(tD);

    cmptMax(tRes.ref(), tD());
    tD.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptMin
(
    meshDirection
    <
        typename meshDirection<Type,MeshType>::cmptType,
        MeshType
    >& res,
    const meshDirection<Type,MeshType>& D
)
{
    cmptMin(res.B(), D.B());
}

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptMin(const meshDirection<Type,MeshType>& D)
{
    typedef typename meshDirection<Type,MeshType>::cmptType cmptType;

    tmp<meshDirection<cmptType,MeshType>> tRes
    (
        new meshDirection<cmptType,MeshType>
        (
            D.fvMsh(),
            D.levelNum(),
            D.directionNum()
        )
    );

    cmptMin(tRes.ref(), D);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptMin(const tmp<meshDirection<Type,MeshType>>& tD)
{
    typedef typename meshDirection<Type,MeshType>::cmptType cmptType;
    tmp<meshDirection<cmptType,MeshType>> tRes =
        reuseDirTmp<cmptType,Type,MeshType>::New(tD);

    cmptMin(tRes.ref(), tD());
    tD.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptAv
(
    meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>& res,
    const meshDirection<Type,MeshType>& D
)
{
    cmptAv(res.B(), D.B());
}

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptAv(const meshDirection<Type,MeshType>& D)
{
    typedef typename meshDirection<Type,MeshType>::cmptType cmptType;

    tmp<meshDirection<cmptType,MeshType>> tRes
    (
        new meshDirection<cmptType,MeshType>
        (
            D.fvMsh(),
            D.levelNum(),
            D.directionNum()
        )
    );

    cmptAv(tRes.ref(), D);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
cmptAv(const tmp<meshDirection<Type,MeshType>>& tD)
{
    typedef typename meshDirection<Type,MeshType>::cmptType cmptType;
    tmp<meshDirection<cmptType,MeshType>> tRes =
        reuseDirTmp<cmptType,Type,MeshType>::New(tD);

    cmptAv(tRes.ref(), tD());
    tD.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptMag
(
    meshDirection<Type,MeshType>& res,
    const meshDirection<Type,MeshType>& D
)
{
    cmptMag(res.B(), D.B());
}

template<class Type, class MeshType>
tmp<meshDirection<Type,MeshType>>
cmptMag(const meshDirection<Type,MeshType>& D)
{
    tmp<meshDirection<Type,MeshType>> tRes
    (
        new meshDirection<Type,MeshType>
        (
            D.fvMsh(),
            D.levelNum(),
            D.directionNum()
        )
    );

    cmptMag(tRes.ref(), D);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshDirection<Type,MeshType>>
cmptMag(const tmp<meshDirection<Type,MeshType>>& tD)
{
    tmp<meshDirection<Type,MeshType>> tRes = New(tD);

    cmptMag(tRes.ref(), tD());
    tD.clear();
    return tRes;
}

template<class Type, class MeshType>
Type max(const meshDirection<Type,MeshType>& D)
{
    Type Max(pTraits<Type>::min);

    // Only for interior cells

    forAllCells(D, i, j, k)
    {
        if (D(i,j,k) > Max)
        {
            Max = D(i,j,k);
        }
    }

    return Max;
}

template<class Type, class MeshType>
Type max(const tmp<meshDirection<Type,MeshType>>& tD)
{
    Type ret(max(tD()));
    tD.clear();
    return ret;
}

template<class Type, class MeshType>
Type min(const meshDirection<Type,MeshType>& D)
{
    Type Min(pTraits<Type>::max);

    // Only for interior cells

    forAllCells(D, i, j, k)
    {
        if (D(i,j,k) < Min)
        {
            Min = D(i,j,k);
        }
    }

    return Min;
}

template<class Type, class MeshType>
Type min(const tmp<meshDirection<Type,MeshType>>& tD)
{
    Type ret(min(tD()));
    tD.clear();
    return ret;
}

template<class Type, class MeshType>
Type sum(const meshDirection<Type,MeshType>& D)
{
    Type Sum(Zero);

    // Only for interior cells

    forAllCells(D, i, j, k)
    {
        Sum += D(i,j,k);
    }

    return Sum;
}

template<class Type, class MeshType>
Type sum(const tmp<meshDirection<Type,MeshType>>& tD)
{
    Type ret(sum(tD()));
    tD.clear();
    return ret;
}

template<class Type, class MeshType>
Type average(const meshDirection<Type,MeshType>& D)
{
    return sum(D)/D.size();
}

template<class Type, class MeshType>
Type average(const tmp<meshDirection<Type,MeshType>>& tD)
{
    Type ret(average(tD()));
    tD.clear();
    return ret;
}

template<class Type, class MeshType>
scalar maxMagSqr(const meshDirection<Type,MeshType>& D)
{
    scalar Max(pTraits<scalar>::min);

    // Only for interior cells

    forAllCells(D, i, j, k)
    {
        if (Foam::magSqr(D(i,j,k)) > Max)
        {
            Max = Foam::magSqr(D(i,j,k));
        }
    }

    return Max;
}

template<class Type, class MeshType>
scalar maxMagSqr(const tmp<meshDirection<Type,MeshType>>& tD)
{
    scalar ret(maxMagSqr(tD()));
    tD.clear();
    return ret;
}

template<class Type, class MeshType>
scalar minMagSqr(const meshDirection<Type,MeshType>& D)
{
    scalar Min(pTraits<scalar>::max);

    // Only for interior cells

    forAllCells(D, i, j, k)
    {
        if (magSqr(D(i,j,k)) < Min)
        {
            Min = Foam::magSqr(D(i,j,k));
        }
    }

    return Min;
}

template<class Type, class MeshType>
scalar minMagSqr(const tmp<meshDirection<Type,MeshType>>& tD)
{
    scalar ret(minMagSqr(tD()));
    tD.clear();
    return ret;
}

template<class Type, class MeshType>
scalar sumProd
(
    const meshDirection<Type,MeshType>& D1,
    const meshDirection<Type,MeshType>& D2
)
{
    scalar SumProd(Zero);

    // Only for interior cells

    forAllCells(D1, i, j, k)
    {
        SumProd += (D1(i,j,k) && D2(i,j,k));
    }

    return SumProd;
}

template<class Type, class MeshType>
scalar sumProd
(
    const tmp<meshDirection<Type,MeshType>>& tD1,
    const meshDirection<Type,MeshType>& D2
)
{
    scalar ret(sumProd(tD1(),D2));
    tD1.clear();
    return ret;
}

template<class Type, class MeshType>
scalar sumProd
(
    const meshDirection<Type,MeshType>& D1,
    const tmp<meshDirection<Type,MeshType>>& tD2
)
{
    scalar ret(sumProd(D1,tD2()));
    tD2.clear();
    return ret;
}

template<class Type, class MeshType>
scalar sumProd
(
    const tmp<meshDirection<Type,MeshType>>& tD1,
    const tmp<meshDirection<Type,MeshType>>& tD2
)
{
    scalar ret(sumProd(tD1(),tD2()));
    tD1.clear();
    tD2.clear();
    return ret;
}

template<class Type, class MeshType>
Type sumCmptProd
(
    const meshDirection<Type,MeshType>& D1,
    const meshDirection<Type,MeshType>& D2
)
{
    Type Sum = Zero;

    // Only for interior cells

    forAllCells(D1, i, j, k)
    {
        Sum += Foam::cmptMultiply(D1(i,j,k),D2(i,j,k));
    }

    return Sum;
}

template<class Type, class MeshType>
Type sumCmptProd
(
    const tmp<meshDirection<Type,MeshType>>& tD1,
    const meshDirection<Type,MeshType>& D2
)
{
    Type ret(sumCmptProd(tD1(),D2));
    tD1.clear();
    return ret;
}

template<class Type, class MeshType>
Type sumCmptProd
(
    const meshDirection<Type,MeshType>& D1,
    const tmp<meshDirection<Type,MeshType>>& tD2
)
{
    Type ret(sumCmptProd(D1,tD2()));
    tD2.clear();
    return ret;
}

template<class Type, class MeshType>
Type sumCmptProd
(
    const tmp<meshDirection<Type,MeshType>>& tD1,
    const tmp<meshDirection<Type,MeshType>>& tD2
)
{
    Type ret(sumCmptProd(tD1(),tD2()));
    tD1.clear();
    tD2.clear();
    return ret;
}

template<class Type, class MeshType>
scalar sumMag(const meshDirection<Type,MeshType>& D)
{
    scalar SumMag(Zero);

    // Only for interior cells

    forAllCells(D, i, j, k)
    {
        SumMag += Foam::mag(D(i,j,k));
    }

    return SumMag;
}

template<class Type, class MeshType>
scalar sumMag(const tmp<meshDirection<Type,MeshType>>& tD)
{
    scalar ret(sumMag(tD()));
    tD.clear();
    return ret;
}

template<class Type, class MeshType>
Type sumCmptMag(const meshDirection<Type,MeshType>& D)
{
    Type Sum(Zero);

    // Only for interior cells

    forAllCells(D, i, j, k)
    {
        Sum += Foam::cmptMag(D(i,j,k));
    }

    return Sum;
}

template<class Type, class MeshType>
Type sumCmptMag(const tmp<meshDirection<Type,MeshType>>& tD)
{
    Type ret(sumCmptMag(tD()));
    tD.clear();
    return ret;
}


#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                        \
                                                                                \
template<class Type, class MeshType>                                            \
ReturnType gFunc                                                                \
(                                                                               \
    const meshDirection<Type,MeshType>& D,                                      \
    const label comm                                                            \
)                                                                               \
{                                                                               \
    ReturnType res(Func(D));                                                    \
    reduce(res, rFunc##Op<ReturnType>(), Pstream::msgType(), comm);             \
    return res;                                                                 \
}                                                                               \
                                                                                \
template<class Type, class MeshType>                                            \
ReturnType gFunc                                                                \
(                                                                               \
    const tmp<meshDirection<Type,MeshType>>& tD,                                \
    const label comm                                                            \
)                                                                               \
{                                                                               \
    ReturnType ret(gFunc(tD()));                                                \
    tD.clear();                                                                 \
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
Type gAverage
(
    const meshDirection<Type,MeshType>& D,
    const label comm
)
{
    label n = D.size();
    Type s = sum(D);
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

template<class Type, class MeshType>
Type gAverage
(
    const tmp<meshDirection<Type,MeshType>>& tD,
    const label comm
)
{
    Type ret(gAverage(tD()));
    tD.clear();
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
