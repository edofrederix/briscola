#include "meshLevelReuseFunctions.H"

#define UNARY_FUNCTION(ReturnType, Type, Func)                                  \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const meshLevel<Type,MeshType>& f                                           \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        Func(res[d], f[d]);                                                     \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func(const meshLevel<Type,MeshType>& f)     \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f.fvMsh(),f.levelNum()));       \
    Func(tRes.ref(), f);                                                        \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>>                                             \
Func(const tmp<meshLevel<Type,MeshType>>& tf)                                   \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type,MeshType>::New(tf);                         \
    Func(tRes.ref(), tf());                                                     \
    tf.clear();                                                                 \
    return tRes;                                                                \
}

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                            \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const meshLevel<Type,MeshType>& f                                           \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d],f[d]);                                                    \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>>                                             \
operator Op(const meshLevel<Type,MeshType>& f)                                  \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f.fvMsh(),f.levelNum()));       \
    OpFunc(tRes.ref(), f);                                                      \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>>                                             \
operator Op(const tmp<meshLevel<Type,MeshType>>& tf)                            \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type,MeshType>::New(tf);                         \
    OpFunc(tRes.ref(), tf());                                                   \
    tf.clear();                                                                 \
    return tRes;                                                                \
}

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                         \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        Func(res[d], f1[d], f2[d]);                                             \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f1.fvMsh(),f1.levelNum()));     \
    Func(tRes.ref(), f1, f2);                                                   \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                   \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type2,MeshType>::New(tf2);                       \
    Func(tRes.ref(), f1, tf2());                                                \
    tf2.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                  \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type1,MeshType>::New(tf1);                       \
    Func(tRes.ref(), tf1(), f2);                                                \
    tf1.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                  \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                   \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmpTmp<ReturnType,Type1,Type1,Type2,MeshType>::New(tf1,tf2);    \
    Func(tRes.ref(), tf1(), tf2());                                             \
    tf1.clear();                                                                \
    tf2.clear();                                                                \
    return tRes;                                                                \
}

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const Type1& s1,                                                            \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        Func(res[d], s1, f2[d]);                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const Type1& s1,                                                            \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f2.fvMsh(),f2.levelNum()));     \
    Func(tRes.ref(), s1, f2);                                                   \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const Type1& s1,                                                            \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                   \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type2,MeshType>::New(tf2);                       \
    Func(tRes.ref(), s1, tf2());                                                \
    tf2.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const List<Type1>& s1,                                                      \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        Func(res[d], s1[d], f2[d]);                                             \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const List<Type1>& s1,                                                      \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f2.fvMsh(),f2.levelNum()));     \
    Func(tRes.ref(), s1, f2);                                                   \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const List<Type1>& s1,                                                      \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                   \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type2,MeshType>::New(tf2);                       \
    Func(tRes.ref(), s1, tf2());                                                \
    tf2.clear();                                                                \
    return tRes;                                                                \
}

#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        Func(res[d], f1[d], s2);                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f1.fvMsh(),f1.levelNum()));     \
    Func(tRes.ref(), f1, s2);                                                   \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                  \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type1,MeshType>::New(tf1);                       \
    Func(tRes.ref(), tf1(), s2);                                                \
    tf1.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const List<Type2>& s2                                                       \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        Func(res[d], f1[d], s2[d]);                                             \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const List<Type2>& s2                                                       \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f1.fvMsh(),f1.levelNum()));     \
    Func(tRes.ref(), f1, s2);                                                   \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                  \
    const List<Type2>& s2                                                       \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type1,MeshType>::New(tf1);                       \
    Func(tRes.ref(), tf1(), s2);                                                \
    tf1.clear();                                                                \
    return tRes;                                                                \
}

#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                     \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                   \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], f1[d], f2[d]);                                           \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f1.fvMsh(),f1.levelNum()));     \
    OpFunc(tRes.ref(), f1, f2);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                   \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type2,MeshType>::New(tf2);                       \
    OpFunc(tRes.ref(), f1, tf2());                                              \
    tf2.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                  \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type1,MeshType>::New(tf1);                       \
    OpFunc(tRes.ref(), tf1(), f2);                                              \
    tf1.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                  \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                   \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmpTmp<ReturnType,Type1,Type1,Type2,MeshType>::New(tf1, tf2);   \
    OpFunc(tRes.ref(), tf1(), tf2());                                           \
    tf1.clear();                                                                \
    tf2.clear();                                                                \
    return tRes;                                                                \
}

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const Type1& s1,                                                            \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], s1, f2[d]);                                              \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const Type1& s1,                                                            \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f2.fvMsh(),f2.levelNum()));     \
    OpFunc(tRes.ref(), s1, f2);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const Type1& s1,                                                            \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                   \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type2,MeshType>::New(tf2);                       \
    OpFunc(tRes.ref(), s1, tf2());                                              \
    tf2.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const List<Type1>& s1,                                                      \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], s1[d], f2[d]);                                           \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const List<Type1>& s1,                                                      \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f2.fvMsh(),f2.levelNum()));     \
    OpFunc(tRes.ref(), s1, f2);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const List<Type1>& s1,                                                      \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                   \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type2,MeshType>::New(tf2);                       \
    OpFunc(tRes.ref(), s1, tf2());                                              \
    tf2.clear();                                                                \
    return tRes;                                                                \
}

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], f1[d], s2);                                              \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f1.fvMsh(),f1.levelNum()));     \
    OpFunc(tRes.ref(), f1, s2);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                  \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type1,MeshType>::New(tf1);                       \
    OpFunc(tRes.ref(), tf1(), s2);                                              \
    tf1.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>& res,                                        \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const List<Type2>& s2                                                       \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], f1[d], s2[d]);                                           \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const List<Type2>& s2                                                       \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>>                                         \
        tRes(new meshLevel<ReturnType,MeshType>(f1.fvMsh(),f1.levelNum()));     \
    OpFunc(tRes.ref(), f1, s2);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                  \
    const List<Type2>& s2                                                       \
)                                                                               \
{                                                                               \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                  \
        reuseLevTmp<ReturnType,Type1,MeshType>::New(tf1);                       \
    OpFunc(tRes.ref(), tf1(), s2);                                              \
    tf1.clear();                                                                \
    return tRes;                                                                \
}

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)               \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                   \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<typename product<Type1, Type2>::type,MeshType>& res,              \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], f1[d], f2[d]);                                           \
}                                                                               \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                   \
operator Op                                                                     \
(                                                                               \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    typedef typename product<Type1, Type2>::type productType;                   \
    tmp<meshLevel<productType,MeshType>>                                        \
        tRes(new meshLevel<productType,MeshType>(f1.fvMsh(),f1.levelNum()));    \
    OpFunc(tRes.ref(), f1, f2);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                   \
operator Op                                                                     \
(                                                                               \
    const meshLevel<Type1,MeshType>& f1,                                        \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                   \
)                                                                               \
{                                                                               \
    typedef typename product<Type1, Type2>::type productType;                   \
    tmp<meshLevel<productType,MeshType>> tRes =                                 \
        reuseLevTmp<productType,Type2,MeshType>::New(tf2);                      \
    OpFunc(tRes.ref(), f1, tf2());                                              \
    tf2.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                   \
operator Op                                                                     \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                  \
    const meshLevel<Type2,MeshType>& f2                                         \
)                                                                               \
{                                                                               \
    typedef typename product<Type1, Type2>::type productType;                   \
    tmp<meshLevel<productType,MeshType>> tRes =                                 \
        reuseLevTmp<productType,Type1,MeshType>::New(tf1);                      \
    OpFunc(tRes.ref(), tf1(), f2);                                              \
    tf1.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshLevel<typename product<Type1, Type2>::type,MeshType>>                   \
operator Op                                                                     \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                  \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                   \
)                                                                               \
{                                                                               \
    typedef typename product<Type1, Type2>::type productType;                   \
    tmp<meshLevel<productType,MeshType>> tRes =                                 \
        reuseLevTmpTmp<productType,Type1,Type1,Type2,MeshType>::New(tf1, tf2);  \
    OpFunc(tRes.ref(), tf1(), tf2());                                           \
    tf1.clear();                                                                \
    tf2.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
/* VectorSpace */                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<typename product<Type, Form>::type,MeshType>& res,                \
    const meshLevel<Type,MeshType>& f1,                                         \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                      \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], f1[d], vs);                                              \
}                                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const meshLevel<Type,MeshType>& f1,                                         \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                      \
)                                                                               \
{                                                                               \
    typedef typename product<Type, Form>::type productType;                     \
    tmp<meshLevel<productType,MeshType>>                                        \
        tRes(new meshLevel<productType,MeshType>(f1.fvMsh(),f1.levelNum()));    \
    OpFunc(tRes.ref(), f1, static_cast<const Form&>(vs));                       \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const tmp<meshLevel<Type,MeshType>>& tf1,                                   \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                      \
)                                                                               \
{                                                                               \
    typedef typename product<Type, Form>::type productType;                     \
    tmp<meshLevel<productType,MeshType>> tRes =                                 \
        reuseLevTmp<productType,Type,MeshType>::New(tf1);                       \
    OpFunc(tRes.ref(), tf1(), static_cast<const Form&>(vs));                    \
    tf1.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<typename product<Form, Type>::type,MeshType>& res,                \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                     \
    const meshLevel<Type,MeshType>& f1                                          \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], vs, f1[d]);                                              \
}                                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                     \
    const meshLevel<Type,MeshType>& f1                                          \
)                                                                               \
{                                                                               \
    typedef typename product<Form, Type>::type productType;                     \
    tmp<meshLevel<productType,MeshType>>                                        \
        tRes(new meshLevel<productType,MeshType>(f1.fvMsh(),f1.levelNum()));    \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), f1);                       \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                     \
    const tmp<meshLevel<Type,MeshType>>& tf1                                    \
)                                                                               \
{                                                                               \
    typedef typename product<Form, Type>::type productType;                     \
    tmp<meshLevel<productType,MeshType>> tRes =                                 \
        reuseLevTmp<productType,Type,MeshType>::New(tf1);                       \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), tf1());                    \
    tf1.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
/* CellSpace */                                                                 \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<typename product<Type, Form>::type,MeshType>& res,                \
    const meshLevel<Type,MeshType>& f1,                                         \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                        \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], f1[d], vs);                                              \
}                                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const meshLevel<Type,MeshType>& f1,                                         \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                        \
)                                                                               \
{                                                                               \
    typedef typename product<Type, Form>::type productType;                     \
    tmp<meshLevel<productType,MeshType>>                                        \
        tRes(new meshLevel<productType,MeshType>(f1.fvMsh(),f1.levelNum()));    \
    OpFunc(tRes.ref(), f1, static_cast<const Form&>(vs));                       \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const tmp<meshLevel<Type,MeshType>>& tf1,                                   \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                        \
)                                                                               \
{                                                                               \
    typedef typename product<Type, Form>::type productType;                     \
    tmp<meshLevel<productType,MeshType>> tRes =                                 \
        reuseLevTmp<productType,Type,MeshType>::New(tf1);                       \
    OpFunc(tRes.ref(), tf1(), static_cast<const Form&>(vs));                    \
    tf1.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<typename product<Form, Type>::type,MeshType>& res,                \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                       \
    const meshLevel<Type,MeshType>& f1                                          \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], vs, f1[d]);                                              \
}                                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                       \
    const meshLevel<Type,MeshType>& f1                                          \
)                                                                               \
{                                                                               \
    typedef typename product<Form, Type>::type productType;                     \
    tmp<meshLevel<productType,MeshType>>                                        \
        tRes(new meshLevel<productType,MeshType>(f1.fvMsh(),f1.levelNum()));    \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), f1);                       \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template                                                                        \
<                                                                               \
    class Type,                                                                 \
    class Form,                                                                 \
    class Cmpt,                                                                 \
    direction nCmpt,                                                            \
    class MeshType                                                              \
>                                                                               \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                       \
    const tmp<meshLevel<Type,MeshType>>& tf1                                    \
)                                                                               \
{                                                                               \
    typedef typename product<Form, Type>::type productType;                     \
    tmp<meshLevel<productType,MeshType>> tRes =                                 \
        reuseLevTmp<productType,Type,MeshType>::New(tf1);                       \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), tf1());                    \
    tf1.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
/* Lists */                                                                     \
                                                                                \
template<class Type, class Form, class MeshType>                                \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<typename product<Type, Form>::type,MeshType>& res,                \
    const meshLevel<Type,MeshType>& f1,                                         \
    const List<Form>& vs                                                        \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], f1[d], vs[d]);                                           \
}                                                                               \
                                                                                \
template<class Type, class Form, class MeshType>                                \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const meshLevel<Type,MeshType>& f1,                                         \
    const List<Form>& vs                                                        \
)                                                                               \
{                                                                               \
    typedef typename product<Type, Form>::type productType;                     \
    tmp<meshLevel<productType,MeshType>>                                        \
        tRes(new meshLevel<productType,MeshType>(f1.fvMsh(),f1.levelNum()));    \
    OpFunc(tRes.ref(), f1, vs);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template<class Type, class Form, class MeshType>                                \
tmp<meshLevel<typename product<Type, Form>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const tmp<meshLevel<Type,MeshType>>& tf1,                                   \
    const List<Form>& vs                                                        \
)                                                                               \
{                                                                               \
    typedef typename product<Type, Form>::type productType;                     \
    tmp<meshLevel<productType,MeshType>> tRes =                                 \
        reuseLevTmp<productType,Type,MeshType>::New(tf1);                       \
    OpFunc(tRes.ref(), tf1(), vs);                                              \
    tf1.clear();                                                                \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template<class Form, class Type, class MeshType>                                \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<typename product<Form, Type>::type,MeshType>& res,                \
    const List<Form>& vs,                                                       \
    const meshLevel<Type,MeshType>& f1                                          \
)                                                                               \
{                                                                               \
    forAll(res, d)                                                              \
        OpFunc(res[d], vs[d], f1[d]);                                           \
}                                                                               \
                                                                                \
template<class Form, class Type, class MeshType>                                \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const List<Form>& vs,                                                       \
    const meshLevel<Type,MeshType>& f1                                          \
)                                                                               \
{                                                                               \
    typedef typename product<Form, Type>::type productType;                     \
    tmp<meshLevel<productType,MeshType>>                                        \
        tRes(new meshLevel<productType,MeshType>(f1.fvMsh(),f1.levelNum()));    \
    OpFunc(tRes.ref(), vs, f1);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template<class Form, class Type, class MeshType>                                \
tmp<meshLevel<typename product<Form, Type>::type,MeshType>>                     \
operator Op                                                                     \
(                                                                               \
    const List<Form>& vs,                                                       \
    const tmp<meshLevel<Type,MeshType>>& tf1                                    \
)                                                                               \
{                                                                               \
    typedef typename product<Form, Type>::type productType;                     \
    tmp<meshLevel<productType,MeshType>> tRes =                                 \
        reuseLevTmp<productType,Type,MeshType>::New(tf1);                       \
    OpFunc(tRes.ref(), vs, tf1());                                              \
    tf1.clear();                                                                \
    return tRes;                                                                \
}
