#include "meshLevelReuseFunctions.H"

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const meshLevel<Type,MeshType>& f                                          \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        Func(res[d], f[d]);                                                    \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func(const meshLevel<Type,MeshType>& f)    \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f.fvMsh(),                                                         \
            f.levelNum()                                                       \
        );                                                                     \
    Func(tRes.ref(), f);                                                       \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>>                                            \
Func(const tmp<meshLevel<Type,MeshType>>& tf)                                  \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type,MeshType>::New(tf);                      \
    Func(tRes.ref(), tf());                                                    \
    if (tf.isTmp()) tf.clear();                                                \
    return tRes;                                                               \
}

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                           \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const meshLevel<Type,MeshType>& f                                          \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d],f[d]);                                                   \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>>                                            \
operator Op(const meshLevel<Type,MeshType>& f)                                 \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f.fvMsh(),                                                         \
            f.levelNum()                                                       \
        );                                                                     \
    OpFunc(tRes.ref(), f);                                                     \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>>                                            \
operator Op(const tmp<meshLevel<Type,MeshType>>& tf)                           \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type,MeshType>::New(tf);                      \
    OpFunc(tRes.ref(), tf());                                                  \
    if (tf.isTmp()) tf.clear();                                                \
    return tRes;                                                               \
}

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        Func(res[d], f1[d], f2[d]);                                            \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
    Func(tRes.ref(), f1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type2,MeshType>::New(tf2);                    \
    Func(tRes.ref(), f1, tf2());                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                 \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type1,MeshType>::New(tf1);                    \
    Func(tRes.ref(), tf1(), f2);                                               \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                 \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmpTmp<ReturnType,Type1,Type1,Type2,MeshType>::              \
        New(tf1,tf2);                                                          \
    Func(tRes.ref(), tf1(), tf2());                                            \
    if (tf1.isTmp()) tf1.clear();                                              \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const Type1& s1,                                                           \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        Func(res[d], s1, f2[d]);                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const Type1& s1,                                                           \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f2.fvMsh(),                                                        \
            f2.levelNum()                                                      \
        );                                                                     \
    Func(tRes.ref(), s1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type2,MeshType>::New(tf2);                    \
    Func(tRes.ref(), s1, tf2());                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const List<Type1>& s1,                                                     \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        Func(res[d], s1[d], f2[d]);                                            \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const List<Type1>& s1,                                                     \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f2.fvMsh(),                                                        \
            f2.levelNum()                                                      \
        );                                                                     \
    Func(tRes.ref(), s1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const List<Type1>& s1,                                                     \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type2,MeshType>::New(tf2);                    \
    Func(tRes.ref(), s1, tf2());                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        Func(res[d], f1[d], s2);                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
    Func(tRes.ref(), f1, s2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                 \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type1,MeshType>::New(tf1);                    \
    Func(tRes.ref(), tf1(), s2);                                               \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        Func(res[d], f1[d], s2[d]);                                            \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
    Func(tRes.ref(), f1, s2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                 \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type1,MeshType>::New(tf1);                    \
    Func(tRes.ref(), tf1(), s2);                                               \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], f1[d], f2[d]);                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
    OpFunc(tRes.ref(), f1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type2,MeshType>::New(tf2);                    \
    OpFunc(tRes.ref(), f1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                 \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type1,MeshType>::New(tf1);                    \
    OpFunc(tRes.ref(), tf1(), f2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                 \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmpTmp<ReturnType,Type1,Type1,Type2,MeshType>::New(tf1, tf2);\
    OpFunc(tRes.ref(), tf1(), tf2());                                          \
    if (tf1.isTmp()) tf1.clear();                                              \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const Type1& s1,                                                           \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], s1, f2[d]);                                             \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const Type1& s1,                                                           \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f2.fvMsh(),                                                        \
            f2.levelNum()                                                      \
        );                                                                     \
    OpFunc(tRes.ref(), s1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type2,MeshType>::New(tf2);                    \
    OpFunc(tRes.ref(), s1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const List<Type1>& s1,                                                     \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], s1[d], f2[d]);                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const List<Type1>& s1,                                                     \
    const meshLevel<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f2.fvMsh(),                                                        \
            f2.levelNum()                                                      \
        );                                                                     \
    OpFunc(tRes.ref(), s1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const List<Type1>& s1,                                                     \
    const tmp<meshLevel<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type2,MeshType>::New(tf2);                    \
    OpFunc(tRes.ref(), s1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], f1[d], s2);                                             \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
    OpFunc(tRes.ref(), f1, s2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                 \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type1,MeshType>::New(tf1);                    \
    OpFunc(tRes.ref(), tf1(), s2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshLevel<ReturnType,MeshType>& res,                                       \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    forAll(res, d)                                                             \
        OpFunc(res[d], f1[d], s2[d]);                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const meshLevel<Type1,MeshType>& f1,                                       \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        meshLevel<ReturnType,MeshType>::New                                    \
        (                                                                      \
            f1.fvMsh(),                                                        \
            f1.levelNum()                                                      \
        );                                                                     \
    OpFunc(tRes.ref(), f1, s2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const tmp<meshLevel<Type1,MeshType>>& tf1,                                 \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    tmp<meshLevel<ReturnType,MeshType>> tRes =                                 \
        reuseLevelTmp<ReturnType,Type1,MeshType>::New(tf1);                    \
    OpFunc(tRes.ref(), tf1(), s2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)
