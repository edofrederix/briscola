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
    tmp<meshField<SCALARPRODTYPE,MeshType>> tRes =
        meshField<SCALARPRODTYPE,MeshType>::New
        (
            "mag("+f.name()+")",
            f.fvMsh()
        );

    tRes->make(f.deep());

    mag(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<SCALARPRODTYPE,MeshType>>
mag(const tmp<meshField<Type,MeshType>>& tf)
{
    tmp<meshField<SCALARPRODTYPE,MeshType>> tRes =
        reuseFieldTmp<SCALARPRODTYPE,Type,MeshType>::New
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
    tmp<meshField<Type,MeshType>> tRes =
        meshField<Type,MeshType>::New
        (
            "cmptMag("+f.name()+")",
            f.fvMsh()
        );

    tRes->make(f.deep());

    cmptMag(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
cmptMag(const tmp<meshField<Type,MeshType>>& tf)
{
    tmp<meshField<Type,MeshType>> tRes =
        reuseFieldTmp<Type,Type,MeshType>::New(tf,"cmptMag("+tf->name()+")");

    cmptMag(tRes.ref(), tf());

    if (tf.isTmp())
        tf.clear();

    return tRes;
}

template<class Type, class MeshType>
void cmptSqr
(
    meshField<Type,MeshType>& res,
    const meshField<Type,MeshType>& f
)
{
    forAll(res, i)
        cmptSqr(res[i],f[i]);
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
cmptSqr(const meshField<Type,MeshType>& f)
{
    tmp<meshField<Type,MeshType>> tRes =
        meshField<Type,MeshType>::New
        (
            "cmptSqr("+f.name()+")",
            f.fvMsh()
        );

    tRes->make(f.deep());

    cmptSqr(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
cmptSqr(const tmp<meshField<Type,MeshType>>& tf)
{
    tmp<meshField<Type,MeshType>> tRes =
        reuseFieldTmp<Type,Type,MeshType>::New(tf,"cmptSqr("+tf->name()+")");

    cmptSqr(tRes.ref(), tf());

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

#define G_UNARY_FUNCTION(ReturnType, gFunc)                                    \
                                                                               \
template<class Type, class MeshType>                                           \
List<ReturnType>                                                               \
gFunc(const meshField<Type,MeshType>& f)                                       \
{                                                                              \
    return gFunc(f[0]);                                                        \
}                                                                              \
                                                                               \
template<class Type, class MeshType>                                           \
List<ReturnType>                                                               \
gFunc(const tmp<meshField<Type,MeshType>>& tf)                                 \
{                                                                              \
    List<ReturnType> ret(gFunc(tf()));                                         \
    if (tf.isTmp()) tf.clear();                                                \
    return ret;                                                                \
}

G_UNARY_FUNCTION(Type, gMax)
G_UNARY_FUNCTION(Type, gMin)
G_UNARY_FUNCTION(Type, gSum)
G_UNARY_FUNCTION(Type, gAverage)

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

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
void OpFunc                                                                    \
(                                                                              \
    meshField<typename product<Type1, Type2>::type, MeshType>& res,            \
    const meshField<Type1,MeshType>& f1,                                       \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], f1[i], f2[i]);                                          \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshField<typename product<Type1, Type2>::type, MeshType>>                 \
operator Op                                                                    \
(                                                                              \
    const meshField<Type1,MeshType>& f1,                                       \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshField<productType,MeshType>> tRes =                                \
        meshField<productType,MeshType>::New                                   \
        (                                                                      \
            "("+f1.name()+#Op+f2.name()+")",                                   \
            f1.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f1.deep() && f2.deep());                                        \
    OpFunc(tRes.ref(), f1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshField<typename product<Type1, Type2>::type, MeshType>>                 \
operator Op                                                                    \
(                                                                              \
    const meshField<Type1,MeshType>& f1,                                       \
    const tmp<meshField<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshField<productType,MeshType>> tRes =                                \
        reuseFieldTmp<productType,Type2,MeshType>::New                         \
        (                                                                      \
            tf2,                                                               \
            "("+f1.name()+#Op+tf2().name()+")"                                 \
        );                                                                     \
    tRes->make(f1.deep() && tf2->deep());                                      \
    OpFunc(tRes.ref(), f1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshField<typename product<Type1, Type2>::type, MeshType>>                 \
operator Op                                                                    \
(                                                                              \
    const tmp<meshField<Type1,MeshType>>& tf1,                                 \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshField<productType,MeshType>> tRes =                                \
        reuseFieldTmp<productType,Type1,MeshType>::New                         \
        (                                                                      \
            tf1,                                                               \
            "("+tf1().name()+#Op+f2.name()+")"                                 \
        );                                                                     \
    tRes->make(tf1->deep() && f2.deep());                                      \
    OpFunc(tRes.ref(), tf1(), f2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshField<typename product<Type1, Type2>::type, MeshType>>                 \
operator Op                                                                    \
(                                                                              \
    const tmp<meshField<Type1,MeshType>>& tf1,                                 \
    const tmp<meshField<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshField<productType,MeshType>> tRes =                                \
        reuseFieldTmpTmp<productType,Type1,Type1,Type2,MeshType>::New          \
        (                                                                      \
            tf1,                                                               \
            tf2,                                                               \
            "("+tf1().name()+#Op+tf2().name()+")"                              \
        );                                                                     \
    tRes->make(tf1->deep() && tf2->deep());                                    \
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
    meshField<typename product<Type, Form>::type,MeshType>& res,               \
    const meshField<Type,MeshType>& f1,                                        \
    const VectorSpace<Form,Cmpt,nCmpt>& v2                                     \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], f1[i], v2);                                             \
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
tmp<meshField<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const meshField<Type,MeshType>& f1,                                        \
    const VectorSpace<Form,Cmpt,nCmpt>& v2                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        meshField<productType,MeshType>::New                                   \
        (                                                                      \
            "("+f1.name()+#Op+Foam::name(v2)+")",                              \
            f1.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f1.deep());                                                     \
    OpFunc(tRes.ref(), f1, v2);                                                \
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
tmp<meshField<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const tmp<meshField<Type,MeshType>>& tf1,                                  \
    const VectorSpace<Form,Cmpt,nCmpt>& v2                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        reuseFieldTmp<productType,Type,MeshType>::New                          \
        (                                                                      \
            tf1,                                                               \
            "("+tf1->name()+#Op+Foam::name(v2)+")"                             \
        );                                                                     \
    tRes->make(tf1->deep());                                                   \
    OpFunc(tRes.ref(), tf1(), v2);                                             \
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
    meshField<typename product<Form,Type>::type,MeshType>& res,                \
    const VectorSpace<Form,Cmpt,nCmpt>& v1,                                    \
    const meshField<Type,MeshType>& f2                                         \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], v1, f2[i]);                                             \
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
tmp<meshField<typename product<Form,Type>::type,MeshType>>                     \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& v1,                                    \
    const meshField<Type,MeshType>& f2                                         \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        meshField<productType,MeshType>::New                                   \
        (                                                                      \
            "("+Foam::name(v1)+#Op+f2.name()+")",                              \
            f2.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f2.deep());                                                     \
    OpFunc(tRes.ref(), v1, f2);                                                \
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
tmp<meshField<typename product<Form,Type>::type,MeshType>>                     \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& v1,                                    \
    const tmp<meshField<Type,MeshType>>& tf2                                   \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        reuseFieldTmp<productType,Type,MeshType>::New                          \
        (                                                                      \
            tf2,                                                               \
            "("+Foam::name(v1)+#Op+tf2->name()+")"                             \
        );                                                                     \
    tRes->make(tf2->deep());                                                   \
    OpFunc(tRes.ref(), v1, tf2);                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
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
    meshField<typename product<Type, Form>::type,MeshType>& res,               \
    const meshField<Type,MeshType>& f1,                                        \
    const CellSpace<Form,Cmpt,nCmpt>& v2                                       \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], f1[i], v2);                                             \
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
tmp<meshField<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const meshField<Type,MeshType>& f1,                                        \
    const CellSpace<Form,Cmpt,nCmpt>& v2                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        meshField<productType,MeshType>::New                                   \
        (                                                                      \
            "("+f1.name()+#Op+Foam::name(v2)+")",                              \
            f1.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f1.deep());                                                     \
    OpFunc(tRes.ref(), f1, v2);                                                \
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
tmp<meshField<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const tmp<meshField<Type,MeshType>>& tf1,                                  \
    const CellSpace<Form,Cmpt,nCmpt>& v2                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        reuseFieldTmp<productType,Type,MeshType>::New                          \
        (                                                                      \
            tf1,                                                               \
            "("+tf1->name()+#Op+Foam::name(v2)+")"                             \
        );                                                                     \
    tRes->make(tf1->deep());                                                   \
    OpFunc(tRes.ref(), tf1(), v2);                                             \
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
    meshField<typename product<Form,Type>::type,MeshType>& res,                \
    const CellSpace<Form,Cmpt,nCmpt>& v1,                                      \
    const meshField<Type,MeshType>& f2                                         \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], v1, f2[i]);                                             \
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
tmp<meshField<typename product<Form,Type>::type,MeshType>>                     \
operator Op                                                                    \
(                                                                              \
    const CellSpace<Form,Cmpt,nCmpt>& v1,                                      \
    const meshField<Type,MeshType>& f2                                         \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        meshField<productType,MeshType>::New                                   \
        (                                                                      \
            "("+Foam::name(v1)+#Op+f2.name()+")",                              \
            f2.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f2.deep());                                                     \
    OpFunc(tRes.ref(), v1, f2);                                                \
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
tmp<meshField<typename product<Form,Type>::type,MeshType>>                     \
operator Op                                                                    \
(                                                                              \
    const CellSpace<Form,Cmpt,nCmpt>& v1,                                      \
    const tmp<meshField<Type,MeshType>>& tf2                                   \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        reuseFieldTmp<productType,Type,MeshType>::New                          \
        (                                                                      \
            tf2,                                                               \
            "("+Foam::name(v1)+#Op+tf2->name()+")"                             \
        );                                                                     \
    tRes->make(tf2->deep());                                                   \
    OpFunc(tRes.ref(), v1, tf2);                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
/* Lists */                                                                    \
                                                                               \
template<class Type, class Form, class MeshType>                               \
void OpFunc                                                                    \
(                                                                              \
    meshField<typename product<Type, Form>::type,MeshType>& res,               \
    const meshField<Type,MeshType>& f1,                                        \
    const List<Form>& v2                                                       \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], f1[i], v2);                                             \
}                                                                              \
                                                                               \
template<class Type, class Form, class MeshType>                               \
tmp<meshField<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const meshField<Type,MeshType>& f1,                                        \
    const List<Form>& v2                                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        meshField<productType,MeshType>::New                                   \
        (                                                                      \
            "("+f1.name()+#Op+Foam::name(v2)+")",                              \
            f1.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f1.deep());                                                     \
    OpFunc(tRes.ref(), f1, v2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type, class Form, class MeshType>                               \
tmp<meshField<typename product<Type, Form>::type,MeshType>>                    \
operator Op                                                                    \
(                                                                              \
    const tmp<meshField<Type,MeshType>>& tf1,                                  \
    const List<Form>& v2                                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        reuseFieldTmp<productType,Type,MeshType>::New                          \
        (                                                                      \
            tf1,                                                               \
            "("+tf1->name()+#Op+Foam::name(v2)+")"                             \
        );                                                                     \
    tRes->make(tf1->deep());                                                   \
    OpFunc(tRes.ref(), tf1(), v2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Form, class Type, class MeshType>                               \
void OpFunc                                                                    \
(                                                                              \
    meshField<typename product<Form,Type>::type,MeshType>& res,                \
    const List<Form>& v1,                                                      \
    const meshField<Type,MeshType>& f2                                         \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], v1, f2[i]);                                             \
}                                                                              \
                                                                               \
template<class Form, class Type, class MeshType>                               \
tmp<meshField<typename product<Form,Type>::type,MeshType>>                     \
operator Op                                                                    \
(                                                                              \
    const List<Form>& v1,                                                      \
    const meshField<Type,MeshType>& f2                                         \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        meshField<productType,MeshType>::New                                   \
        (                                                                      \
            "("+Foam::name(v1)+#Op+f2.name()+")",                              \
            f2.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f2.deep());                                                     \
    OpFunc(tRes.ref(), v1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Form, class Type, class MeshType>                               \
tmp<meshField<typename product<Form,Type>::type,MeshType>>                     \
operator Op                                                                    \
(                                                                              \
    const List<Form>& v1,                                                      \
    const tmp<meshField<Type,MeshType>>& tf2                                   \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshField<productType,MeshType>> tRes =                                \
        reuseFieldTmp<productType,Type,MeshType>::New                          \
        (                                                                      \
            tf2,                                                               \
            "("+Foam::name(v1)+#Op+tf2->name()+")"                             \
        );                                                                     \
    tRes->make(tf2->deep());                                                   \
    OpFunc(tRes.ref(), v1, tf2);                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
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
