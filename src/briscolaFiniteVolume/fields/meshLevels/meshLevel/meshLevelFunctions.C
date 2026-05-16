#include "PstreamReduceOps.H"
#include "meshLevelReuseFunctions.H"

#define TEMPLATE template<class Type, class MeshType>
#include "meshLevelFunctionsM.C"

// Scalar return type must be deduced because of cell space
#define SCALARPRODTYPE typename scalarProduct<Type,Type>::type

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void mag(meshLevel<SCALARPRODTYPE,MeshType>& res, const meshLevel<Type,MeshType>& f)
{
    forAll(res, d)
        mag(res[d], f[d]);
}

template<class Type, class MeshType>
tmp<meshLevel<SCALARPRODTYPE,MeshType>> mag(const meshLevel<Type,MeshType>& f)
{
    tmp<meshLevel<SCALARPRODTYPE,MeshType>> tRes =
        meshLevel<SCALARPRODTYPE,MeshType>::New
        (
            f.fvMsh(),
            f.levelNum()
        );

    mag(tRes.ref(), f);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshLevel<SCALARPRODTYPE,MeshType>> mag(const tmp<meshLevel<Type,MeshType>>& tf)
{
    tmp<meshLevel<SCALARPRODTYPE,MeshType>> tRes =
        reuseLevelTmp<SCALARPRODTYPE,Type,MeshType>::New(tf);

    mag(tRes.ref(), tf());
    if (tf.isTmp())
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
    tmp<meshLevel<Type,MeshType>> tRes =
        meshLevel<Type,MeshType>::New
        (
            f.fvMsh(),
            f.levelNum()
        );

    cmptMag(tRes.ref(), f);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptMag(const tmp<meshLevel<Type,MeshType>>& tf)
{
    tmp<meshLevel<Type,MeshType>> tRes =
        reuseLevelTmp<Type,Type,MeshType>::New(tf);

    cmptMag(tRes.ref(), tf());
    if (tf.isTmp())
        tf.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptSqr(meshLevel<Type,MeshType>& res, const meshLevel<Type,MeshType>& f)
{
    forAll(res, d)
        cmptSqr(res[d],f[d]);
}

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptSqr(const meshLevel<Type,MeshType>& f)
{
    tmp<meshLevel<Type,MeshType>> tRes =
        meshLevel<Type,MeshType>::New
        (
            f.fvMsh(),
            f.levelNum()
        );

    cmptSqr(tRes.ref(), f);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshLevel<Type,MeshType>> cmptSqr(const tmp<meshLevel<Type,MeshType>>& tf)
{
    tmp<meshLevel<Type,MeshType>> tRes =
        reuseLevelTmp<Type,Type,MeshType>::New(tf);

    cmptSqr(tRes.ref(), tf());
    if (tf.isTmp())
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
    if (tf.isTmp())
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
    if (tf.isTmp())
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
    if (tf.isTmp())
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
    if (tf.isTmp())
        tf.clear();
    return ret;
}

#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                       \
                                                                               \
template<class Type, class MeshType>                                           \
List<ReturnType> gFunc(const meshLevel<Type,MeshType>& f)                      \
{                                                                              \
    List<ReturnType> res(Func(f));                                             \
                                                                               \
    /* We need to communicate per direction, because not all functions are */  \
    /* defined for lists */                                                    \
                                                                               \
    forAll(res, i)                                                             \
        reduce(res[i], rFunc##Op<ReturnType>(), Pstream::msgType());           \
                                                                               \
    return res;                                                                \
}                                                                              \
                                                                               \
template<class Type, class MeshType>                                           \
List<ReturnType> gFunc(const tmp<meshLevel<Type,MeshType>>& tf)                \
{                                                                              \
    List<ReturnType> ret(gFunc(tf()));                                         \
    if (tf.isTmp()) tf.clear();                                                \
    return ret;                                                                \
}

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)

#undef G_UNARY_FUNCTION

template<class Type, class MeshType>
List<Type> gAverage(const meshLevel<Type,MeshType>& f)
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

        sumReduce(Sum[d], n, Pstream::msgType(), f.lvl().comms());

        if (n > 0)
        {
            Avrg[d] = Sum[d]/n;
        }
    }

    return Avrg;
}

template<class Type, class MeshType>
List<Type> gAverage(const tmp<meshLevel<Type,MeshType>>& tf)
{
    List<Type> ret(gAverage(tf()));
    if (tf.isTmp())
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


#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<typename product<Type1, Type2>::type,MeshType>& res,             \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], f1[d], f2[d]);                                          \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                  \
operator Op                                                                    \
(                                                                              \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        meshLevel<productType,MeshType>::New                                   \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
    OpFunc(tRes.ref(), f1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                  \
operator Op                                                                    \
(                                                                              \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        reuseLevelTmp<productType,Type2,MeshType>::New(tf2);                   \
    OpFunc(tRes.ref(), f1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                  \
operator Op                                                                    \
(                                                                              \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                 \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        reuseLevelTmp<productType,Type1,MeshType>::New(tf1);                   \
    OpFunc(tRes.ref(), tf1(), f2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                  \
operator Op                                                                    \
(                                                                              \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                 \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        reuseLevelTmpTmp<productType,Type1,Type1,Type2,MeshType>::             \
        New(tf1, tf2);                                                         \
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
    class MeshType                                                             \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<typename product<Type, Form>::type,MeshType>& res,               \
    const meshLevel<Type,MeshType>& f1,                                        \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], f1[d], vs);                                             \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const meshLevel<Type,MeshType>& f1,                                        \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        meshLevel<productType,MeshType>::New                                   \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
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
    class MeshType                                                             \
>                                                                              \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const tmp<meshLevel<Type,MeshType>>& tf1,                                  \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        reuseLevelTmp<productType,Type,MeshType>::New(tf1);                    \
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
    class MeshType                                                             \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<typename product<Form, Type>::type,MeshType>& res,               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const meshLevel<Type,MeshType>& f1                                         \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], vs, f1[d]);                                             \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const meshLevel<Type,MeshType>& f1                                         \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        meshLevel<productType,MeshType>::New                                   \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
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
    class MeshType                                                             \
>                                                                              \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<meshLevel<Type,MeshType>>& tf1                                   \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        reuseLevelTmp<productType,Type,MeshType>::New(tf1);                    \
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
    class MeshType                                                             \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<typename product<Type, Form>::type,MeshType>& res,               \
    const meshLevel<Type,MeshType>& f1,                                        \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], f1[d], vs);                                             \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const meshLevel<Type,MeshType>& f1,                                        \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        meshLevel<productType,MeshType>::New                                   \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
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
    class MeshType                                                             \
>                                                                              \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const tmp<meshLevel<Type,MeshType>>& tf1,                                  \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        reuseLevelTmp<productType,Type,MeshType>::New(tf1);                    \
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
    class MeshType                                                             \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<typename product<Form, Type>::type,MeshType>& res,               \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const meshLevel<Type,MeshType>& f1                                         \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], vs, f1[d]);                                             \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const meshLevel<Type,MeshType>& f1                                         \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        meshLevel<productType,MeshType>::New                                   \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
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
    class MeshType                                                             \
>                                                                              \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const tmp<meshLevel<Type,MeshType>>& tf1                                   \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        reuseLevelTmp<productType,Type,MeshType>::New(tf1);                    \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), tf1());                   \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
/* Lists */                                                                    \
                                                                               \
template<class Type, class Form, class MeshType>                               \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<typename product<Type, Form>::type,MeshType>& res,               \
    const meshLevel<Type,MeshType>& f1,                                        \
    const List<Form>& vs                                                       \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], f1[d], vs[d]);                                          \
}                                                                              \
                                                                               \
template<class Type, class Form, class MeshType>                               \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const meshLevel<Type,MeshType>& f1,                                        \
    const List<Form>& vs                                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        meshLevel<productType,MeshType>::New                                   \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
    OpFunc(tRes.ref(), f1, vs);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type, class Form, class MeshType>                               \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const tmp<meshLevel<Type,MeshType>>& tf1,                                  \
    const List<Form>& vs                                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        reuseLevelTmp<productType,Type,MeshType>::New(tf1);                    \
    OpFunc(tRes.ref(), tf1(), vs);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Form, class Type, class MeshType>                               \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<typename product<Form, Type>::type,MeshType>& res,               \
    const List<Form>& vs,                                                      \
    const meshLevel<Type,MeshType>& f1                                         \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], vs[d], f1[d]);                                          \
}                                                                              \
                                                                               \
template<class Form, class Type, class MeshType>                               \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const List<Form>& vs,                                                      \
    const meshLevel<Type,MeshType>& f1                                         \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        meshLevel<productType,MeshType>::New                                   \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
    OpFunc(tRes.ref(), vs, f1);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Form, class Type, class MeshType>                               \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const List<Form>& vs,                                                      \
    const tmp<meshLevel<Type,MeshType>>& tf1                                   \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshLevel<productType,MeshType>> tRes =                                \
        reuseLevelTmp<productType,Type,MeshType>::New(tf1);                    \
    OpFunc(tRes.ref(), vs, tf1());                                             \
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

}

#undef SCALARPRODTYPE

#include "undefBlockFunctionsM.H"
