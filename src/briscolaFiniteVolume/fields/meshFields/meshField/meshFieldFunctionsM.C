#include "meshFieldReuseFunctions.H"

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshField<ReturnType,MeshType>& res,                                       \
    const meshField<Type,MeshType>& f                                          \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        Func(res[i], f[i]);                                                    \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func(const meshField<Type,MeshType>& f)    \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            #Func "("+f.name()+")",                                            \
            f.fvMsh()                                                          \
        );                                                                     \
    tRes->make(f.deep());                                                      \
    Func(tRes.ref(), f);                                                       \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const tmp<meshField<Type,MeshType>>& tf                                    \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type,MeshType>::New                           \
        (                                                                      \
            tf,                                                                \
            #Func "("+tf().name()+")"                                          \
        );                                                                     \
    Func(tRes.ref(), tf());                                                    \
    if (tf.isTmp()) tf.clear();                                                \
    return tRes;                                                               \
}

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                           \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshField<ReturnType,MeshType>& res,                                       \
    const meshField<Type,MeshType>& f                                          \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], f[i]);                                                  \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const meshField<Type,MeshType>& f                                          \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            #Op "("+f.name()+")",                                              \
            f.fvMsh()                                                          \
        );                                                                     \
    tRes->make(f.deep());                                                      \
    OpFunc(tRes.ref(), f);                                                     \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const tmp<meshField<Type,MeshType>>& tf                                    \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type,MeshType>::New                           \
        (                                                                      \
            tf,                                                                \
            #Op "("+tf().name()+")"                                            \
        );                                                                     \
    OpFunc(tRes.ref(), tf());                                                  \
    if (tf.isTmp()) tf.clear();                                                \
    return tRes;                                                               \
}

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshField<ReturnType,MeshType>& res,                                       \
    const meshField<Type1,MeshType>& f1,                                       \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        Func(res[i], f1[i], f2[i]);                                            \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const meshField<Type1,MeshType>& f1,                                       \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            #Func "("+f1.name()+","+f2.name()+")",                             \
            f1.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f1.deep() && f2.deep());                                        \
    Func(tRes.ref(), f1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const meshField<Type1,MeshType>& f1,                                       \
    const tmp<meshField<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type2,MeshType>::New                          \
        (                                                                      \
            tf2,                                                               \
            #Func "("+f1.name()+","+tf2().name()+")"                           \
        );                                                                     \
    tRes->make(f1.deep() && tf2->deep());                                      \
    Func(tRes.ref(), f1, tf2());                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const tmp<meshField<Type1,MeshType>>& tf1,                                 \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type1,MeshType>::New                          \
        (                                                                      \
            tf1,                                                               \
            #Func "("+tf1().name()+","+f2.name()+")"                           \
        );                                                                     \
    tRes->make(tf1->deep() && f2.deep());                                      \
    Func(tRes.ref(), tf1(), f2);                                               \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const tmp<meshField<Type1,MeshType>>& tf1,                                 \
    const tmp<meshField<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmpTmp<ReturnType,Type1,Type1,Type2,MeshType>::New           \
        (                                                                      \
            tf1,                                                               \
            tf2,                                                               \
            #Func "("+tf1().name()+","+tf2().name()+")"                        \
        );                                                                     \
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
    meshField<ReturnType,MeshType>& res,                                       \
    const Type1& s1,                                                           \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        Func(res[i], s1, f2[i]);                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const Type1& s1,                                                           \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            #Func "("+Foam::name(s1)+","+f2.name()+")",                        \
            f2.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f2.deep());                                                     \
    Func(tRes.ref(), s1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<meshField<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type2,MeshType>::New                          \
        (                                                                      \
            tf2,                                                               \
            #Func "("+Foam::name(s1)+","+tf2().name()+")"                      \
        );                                                                     \
    Func(tRes.ref(), s1, tf2());                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshField<ReturnType,MeshType>& res,                                       \
    const List<Type1>& s1,                                                     \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        Func(res[i], s1, f2[i]);                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const List<Type1>& s1,                                                     \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            #Func "("+Foam::name(s1)+","+f2.name()+")",                        \
            f2.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f2.deep());                                                     \
    Func(tRes.ref(), s1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const List<Type1>& s1,                                                     \
    const tmp<meshField<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type2,MeshType>::New                          \
        (                                                                      \
            tf2,                                                               \
            #Func "("+Foam::name(s1)+","+tf2().name()+")"                      \
        );                                                                     \
    Func(tRes.ref(), s1, tf2());                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshField<ReturnType,MeshType>& res,                                       \
    const meshField<Type1,MeshType>& f1,                                       \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        Func(res[i], f1[i], s2);                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const meshField<Type1,MeshType>& f1,                                       \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            #Func "("+f1.name()+","+Foam::name(s2)+")",                        \
            f1.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f1.deep());                                                     \
    Func(tRes.ref(), f1, s2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const tmp<meshField<Type1,MeshType>>& tf1,                                 \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type1,MeshType>::New                          \
        (                                                                      \
            tf1,                                                               \
            #Func "("+tf1().name()+","+Foam::name(s2)+")"                      \
        );                                                                     \
    Func(tRes.ref(), tf1(), s2);                                               \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshField<ReturnType,MeshType>& res,                                       \
    const meshField<Type1,MeshType>& f1,                                       \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        Func(res[i], f1[i], s2);                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const meshField<Type1,MeshType>& f1,                                       \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            #Func "("+f1.name()+","+Foam::name(s2)+")",                        \
            f1.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f1.deep());                                                     \
    Func(tRes.ref(), f1, s2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> Func                                       \
(                                                                              \
    const tmp<meshField<Type1,MeshType>>& tf1,                                 \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type1,MeshType>::New                          \
        (                                                                      \
            tf1,                                                               \
            #Func "("+tf1().name()+","+Foam::name(s2)+")"                      \
        );                                                                     \
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
    meshField<ReturnType,MeshType>& res,                                       \
    const meshField<Type1,MeshType>& f1,                                       \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], f1[i], f2[i]);                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const meshField<Type1,MeshType>& f1,                                       \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            "("+f1.name()+#Op+f2.name()+")",                                   \
            f1.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f1.deep() && f2.deep());                                        \
    OpFunc(tRes.ref(), f1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const meshField<Type1,MeshType>& f1,                                       \
    const tmp<meshField<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type2,MeshType>::New                          \
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
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const tmp<meshField<Type1,MeshType>>& tf1,                                 \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type1,MeshType>::New                          \
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
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const tmp<meshField<Type1,MeshType>>& tf1,                                 \
    const tmp<meshField<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmpTmp<ReturnType,Type1,Type1,Type2,MeshType>::New           \
        (                                                                      \
            tf1,                                                               \
            tf2,                                                               \
            "("+tf1().name()+#Op+tf2().name()+")"                              \
        );                                                                     \
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
    meshField<ReturnType,MeshType>& res,                                       \
    const Type1& s1,                                                           \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], s1, f2[i]);                                             \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const Type1& s1,                                                           \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            "("+Foam::name(s1)+#Op+f2.name()+")",                              \
            f2.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f2.deep());                                                     \
    OpFunc(tRes.ref(), s1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<meshField<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type2,MeshType>::New                          \
        (                                                                      \
            tf2,                                                               \
            "("+Foam::name(s1)+#Op+tf2().name()+")"                            \
        );                                                                     \
    OpFunc(tRes.ref(), s1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshField<ReturnType,MeshType>& res,                                       \
    const List<Type1>& s1,                                                     \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], s1, f2[i]);                                             \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const List<Type1>& s1,                                                     \
    const meshField<Type2,MeshType>& f2                                        \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            "("+Foam::name(s1)+#Op+f2.name()+")",                              \
            f2.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f2.deep());                                                     \
    OpFunc(tRes.ref(), s1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const List<Type1>& s1,                                                     \
    const tmp<meshField<Type2,MeshType>>& tf2                                  \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type2,MeshType>::New                          \
        (                                                                      \
            tf2,                                                               \
            "("+Foam::name(s1)+#Op+tf2().name()+")"                            \
        );                                                                     \
    OpFunc(tRes.ref(), s1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshField<ReturnType,MeshType>& res,                                       \
    const meshField<Type1,MeshType>& f1,                                       \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], f1[i], s2);                                             \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const meshField<Type1,MeshType>& f1,                                       \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            "("+f1.name()+#Op+Foam::name(s2)+")",                              \
            f1.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f1.deep());                                                     \
    OpFunc(tRes.ref(), f1, s2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const tmp<meshField<Type1,MeshType>>& tf1,                                 \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type1,MeshType>::New                          \
        (                                                                      \
            tf1,                                                               \
            "("+tf1().name()+#Op+Foam::name(s2)+")"                            \
        );                                                                     \
    OpFunc(tRes.ref(), tf1(), s2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshField<ReturnType,MeshType>& res,                                       \
    const meshField<Type1,MeshType>& f1,                                       \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
        OpFunc(res[i], f1[i], s2);                                             \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const meshField<Type1,MeshType>& f1,                                       \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        meshField<ReturnType,MeshType>::New                                    \
        (                                                                      \
            "("+f1.name()+#Op+Foam::name(s2)+")",                              \
            f1.fvMsh()                                                         \
        );                                                                     \
    tRes->make(f1.deep());                                                     \
    OpFunc(tRes.ref(), f1, s2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshField<ReturnType,MeshType>> operator Op                                \
(                                                                              \
    const tmp<meshField<Type1,MeshType>>& tf1,                                 \
    const List<Type2>& s2                                                      \
)                                                                              \
{                                                                              \
    tmp<meshField<ReturnType,MeshType>> tRes =                                 \
        reuseFieldTmp<ReturnType,Type1,MeshType>::New                          \
        (                                                                      \
            tf1,                                                               \
            "("+tf1().name()+#Op+Foam::name(s2)+")"                            \
        );                                                                     \
    OpFunc(tRes.ref(), tf1(), s2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)

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
    OpFunc(tRes.ref(), v1, tf2);                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}
