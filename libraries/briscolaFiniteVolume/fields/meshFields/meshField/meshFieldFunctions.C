#include "PstreamReduceOps.H"
#include "meshFieldReuseFunctions.H"

#define TEMPLATE template<class Type, class MeshType>
#include "meshFieldFunctionsM.C"

#define SCALARPRODTYPE typename scalarProduct<Type,Type>::type

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void mag
(
    meshField<SCALARPRODTYPE,MeshType>& res,
    const meshField<Type,MeshType>& f
)
{
    forAll(res, i)
        mag(res[i],f[i]);
}

template<class Type, class MeshType>
tmp<meshField<SCALARPRODTYPE,MeshType>>
mag(const meshField<Type,MeshType>& f)
{
    tmp<meshField<SCALARPRODTYPE,MeshType>> tRes
    (
        new meshField<SCALARPRODTYPE,MeshType>
        (
            "mag("+f.name()+")",
            f.fvMsh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false,
            f.deep()
        )
    );

    mag(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<SCALARPRODTYPE,MeshType>>
mag(const tmp<meshField<Type,MeshType>>& tf)
{
    tmp<meshField<SCALARPRODTYPE,MeshType>> tRes =
        reuseFieTmp<SCALARPRODTYPE,Type,MeshType>::New
        (
            tf,
            "mag("+tf->name()+")"
        );

    mag(tRes.ref(), tf());

    if (tf.isTmp())
        tf.clear();

    return tRes;
}

template<class Type, class MeshType>
void cmptMax
(
    meshField<typename meshField<Type,MeshType>::cmptType,MeshType>& res,
    const meshField<Type,MeshType>& f
)
{
    forAll(res, i)
        cmptMax(res[i],f[i]);
}

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptMax(const meshField<Type,MeshType>& f)
{
    tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>> tRes
    (
        new meshField<typename meshField<Type,MeshType>::cmptType,MeshType>
        (
            "cmptMax("+f.name()+")",
            f.fvMsh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false,
            f.deep()
        )
    );

    cmptMax(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptMax(const tmp<meshField<Type,MeshType>>& tf)
{
    typedef typename meshField<Type,MeshType>::cmptType cmptType;
    tmp<meshField<cmptType,MeshType>> tRes =
        reuseFieTmp<cmptType,Type,MeshType>::New(tf,"cmptMax("+tf->name()+")");

    cmptMax(tRes.ref(), tf());

    if (tf.isTmp())
        tf.clear();

    return tRes;
}

template<class Type, class MeshType>
void cmptMin
(
    meshField<typename meshField<Type,MeshType>::cmptType,MeshType>& res,
    const meshField<Type,MeshType>& f
)
{
    forAll(res, i)
        cmptMin(res[i],f[i]);
}

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptMin(const meshField<Type,MeshType>& f)
{
    tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>> tRes
    (
        new meshField<typename meshField<Type,MeshType>::cmptType,MeshType>
        (
            "cmptMin("+f.name()+")",
            f.fvMsh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false,
            f.deep()
        )
    );

    cmptMin(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptMin(const tmp<meshField<Type,MeshType>>& tf)
{
    typedef typename meshField<Type,MeshType>::cmptType cmptType;
    tmp<meshField<cmptType,MeshType>> tRes =
        reuseFieTmp<cmptType,Type,MeshType>::New(tf,"cmptMin("+tf->name()+")");

    cmptMin(tRes.ref(), tf());

    if (tf.isTmp())
        tf.clear();

    return tRes;
}

template<class Type, class MeshType>
void cmptAv
(
    meshField<typename meshField<Type,MeshType>::cmptType,MeshType>& res,
    const meshField<Type,MeshType>& f
)
{
    forAll(res, i)
        cmptAv(res[i],f[i]);
}

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptAv(const meshField<Type,MeshType>& f)
{
    tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>> tRes
    (
        new meshField<typename meshField<Type,MeshType>::cmptType,MeshType>
        (
            "cmptAv("+f.name()+")",
            f.fvMsh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false,
            f.deep()
        )
    );

    cmptAv(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<typename meshField<Type,MeshType>::cmptType,MeshType>>
cmptAv(const tmp<meshField<Type,MeshType>>& tf)
{
    typedef typename meshField<Type,MeshType>::cmptType cmptType;
    tmp<meshField<cmptType,MeshType>> tRes =
        reuseFieTmp<cmptType,Type,MeshType>::New(tf,"cmptAv("+tf->name()+")");

    cmptAv(tRes.ref(), tf());

    if (tf.isTmp())
        tf.clear();

    return tRes;
}

template<class Type, class MeshType>
void cmptMag
(
    meshField<Type,MeshType>& res,
    const meshField<Type,MeshType>& f
)
{
    forAll(res, i)
        cmptMag(res[i],f[i]);
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
cmptMag(const meshField<Type,MeshType>& f)
{
    tmp<meshField<Type,MeshType>> tRes
    (
        new meshField<Type,MeshType>
        (
            "cmptMag("+f.name()+")",
            f.fvMsh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false,
            f.deep()
        )
    );

    cmptMag(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
cmptMag(const tmp<meshField<Type,MeshType>>& tf)
{
    tmp<meshField<Type,MeshType>> tRes = New(tf,"cmptMag("+tf->name()+")");

    cmptMag(tRes.ref(), tf());

    if (tf.isTmp())
        tf.clear();

    return tRes;
}

template<class Type, class MeshType>
List<Type>
max(const meshField<Type,MeshType>& f)
{
    return max(f[0]);
}

template<class Type, class MeshType>
List<Type>
max(const tmp<meshField<Type,MeshType>>& tf)
{
    List<Type> ret(max(tf()[0]));
    if (tf.isTmp())
        tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type>
min(const meshField<Type,MeshType>& f)
{
    return min(f[0]);
}

template<class Type, class MeshType>
List<Type>
min(const tmp<meshField<Type,MeshType>>& tf)
{
    List<Type> ret(min(tf()[0]));
    if (tf.isTmp())
        tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type>
sum(const meshField<Type,MeshType>& f)
{
    return sum(f[0]);
}

template<class Type, class MeshType>
List<Type>
sum(const tmp<meshField<Type,MeshType>>& tf)
{
    List<Type> ret(sum(tf()[0]));
    if (tf.isTmp())
        tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type>
average(const meshField<Type,MeshType>& f)
{
    return average(f[0]);
}

template<class Type, class MeshType>
List<Type>
average(const tmp<meshField<Type,MeshType>>& tf)
{
    List<Type> ret(average(tf()[0]));
    if (tf.isTmp())
        tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<SCALARPRODTYPE>
maxMagSqr(const meshField<Type,MeshType>& f)
{
    return maxMagSqr(f[0]);
}

template<class Type, class MeshType>
List<SCALARPRODTYPE>
maxMagSqr(const tmp<meshField<Type,MeshType>>& tf)
{
    List<SCALARPRODTYPE> ret(maxMagSqr(tf()[0]));
    if (tf.isTmp())
        tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<SCALARPRODTYPE>
minMagSqr(const meshField<Type,MeshType>& f)
{
    return minMagSqr(f[0]);
}

template<class Type, class MeshType>
List<SCALARPRODTYPE>
minMagSqr(const tmp<meshField<Type,MeshType>>& tf)
{
    List<SCALARPRODTYPE> ret(minMagSqr(tf()[0]));
    if (tf.isTmp())
        tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<SCALARPRODTYPE>
sumProd
(
    const meshField<Type,MeshType>& f1,
    const meshField<Type,MeshType>& f2
)
{
    return sumProd(f1[0], f2[0]);
}

template<class Type, class MeshType>
List<SCALARPRODTYPE>
sumProd
(
    const tmp<meshField<Type,MeshType>>& tf1,
    const meshField<Type,MeshType>& f2
)
{
    List<SCALARPRODTYPE> ret(sumProd(tf1()[0], f2[0]));
    if (tf1.isTmp())
        tf1.clear();
    return ret;
}

template<class Type, class MeshType>
List<SCALARPRODTYPE>
sumProd
(
    const meshField<Type,MeshType>& f1,
    const tmp<meshField<Type,MeshType>>& tf2
)
{
    List<SCALARPRODTYPE> ret(sumProd(f1[0], tf2()[0]));
    if (tf2.isTmp())
        tf2.clear();
    return ret;
}

template<class Type, class MeshType>
List<SCALARPRODTYPE>
sumProd
(
    const tmp<meshField<Type,MeshType>>& tf1,
    const tmp<meshField<Type,MeshType>>& tf2
)
{
    List<SCALARPRODTYPE> ret(sumProd(tf1()[0], tf2()[0]));
    if (tf1.isTmp())
        tf1.clear();
    if (tf2.isTmp())
        tf2.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const meshField<Type,MeshType>& f1,
    const meshField<Type,MeshType>& f2
)
{
    return sumCmptProd(f1[0], f2[0]);
}

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const tmp<meshField<Type,MeshType>>& tf1,
    const meshField<Type,MeshType>& f2
)
{
    List<Type> ret(sumCmptProd(tf1()[0], f2[0]));
    if (tf1.isTmp())
        tf1.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const meshField<Type,MeshType>& f1,
    const tmp<meshField<Type,MeshType>>& tf2
)
{
    List<Type> ret(sumCmptProd(f1[0], tf2()[0]));
    if (tf2.isTmp())
        tf2.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type> sumCmptProd
(
    const tmp<meshField<Type,MeshType>>& tf1,
    const tmp<meshField<Type,MeshType>>& tf2
)
{
    List<Type> ret(sumCmptProd(tf1()[0], tf2()[0]));
    if (tf1.isTmp())
        tf1.clear();
    if (tf2.isTmp())
        tf2.clear();
    return ret;
}

template<class Type, class MeshType>
List<SCALARPRODTYPE>
sumMag(const meshField<Type,MeshType>& f)
{
    return sumMag(f[0]);
}

template<class Type, class MeshType>
List<SCALARPRODTYPE>
sumMag(const tmp<meshField<Type,MeshType>>& tf)
{
    List<SCALARPRODTYPE> ret(sumMag(tf()[0]));
    if (tf.isTmp())
        tf.clear();
    return ret;
}

template<class Type, class MeshType>
List<Type>
sumCmptMag(const meshField<Type,MeshType>& f)
{
    return sumCmptMag(f[0]);
}

template<class Type, class MeshType>
List<Type>
sumCmptMag(const tmp<meshField<Type,MeshType>>& tf)
{
    List<Type> ret(sumCmptMag(tf()[0]));
    if (tf.isTmp())
        tf.clear();
    return ret;
}

#define G_UNARY_FUNCTION(ReturnType, gFunc)                                    \
                                                                               \
template<class Type, class MeshType>                                           \
List<ReturnType>                                                               \
gFunc(const meshField<Type,MeshType>& f, const label comm)                     \
{                                                                              \
    return gFunc(f[0], comm);                                                  \
}                                                                              \
                                                                               \
template<class Type, class MeshType>                                           \
List<ReturnType>                                                               \
gFunc(const tmp<meshField<Type,MeshType>>& tf, const label comm)               \
{                                                                              \
    List<ReturnType> ret(gFunc(tf->operator[](0),comm));                       \
    if (tf.isTmp()) tf.clear();                                                \
    return ret;                                                                \
}

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
