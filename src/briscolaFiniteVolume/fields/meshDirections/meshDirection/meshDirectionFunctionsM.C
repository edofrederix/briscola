#include "meshDirectionReuseFunctions.H"

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshDirection<ReturnType,MeshType>& res,                                   \
    const meshDirection<Type,MeshType>& D                                      \
)                                                                              \
{                                                                              \
    Func(res.B(), D.B());                                                      \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>>                                        \
Func(const meshDirection<Type,MeshType>& D)                                    \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        meshDirection<ReturnType,MeshType>::New                                \
        (                                                                      \
            D.fvMsh(),                                                         \
            D.levelNum(),                                                      \
            D.directionNum()                                                   \
        );                                                                     \
    Func(tRes.ref(), D);                                                       \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>>                                        \
Func(const tmp<meshDirection<Type,MeshType>>& tD)                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmp<ReturnType,Type,MeshType>::New(tD);                  \
    Func(tRes.ref(), tD());                                                    \
    if (tD.isTmp()) tD.clear();                                                \
    return tRes;                                                               \
}

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                           \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshDirection<ReturnType,MeshType>& res,                                   \
    const meshDirection<Type,MeshType>& D                                      \
)                                                                              \
{                                                                              \
    OpFunc(res.B(), D.B());                                                    \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>>                                        \
operator Op(const meshDirection<Type,MeshType>& D)                             \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        meshDirection<ReturnType,MeshType>::New                                \
        (                                                                      \
            D.fvMsh(),                                                         \
            D.levelNum(),                                                      \
            D.directionNum()                                                   \
        );                                                                     \
    OpFunc(tRes.ref(), D);                                                     \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>>                                        \
operator Op(const tmp<meshDirection<Type,MeshType>>& tD)                       \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmp<ReturnType,Type,MeshType>::New(tD);                  \
    OpFunc(tRes.ref(), tD());                                                  \
    if (tD.isTmp()) tD.clear();                                                \
    return tRes;                                                               \
}

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshDirection<ReturnType,MeshType>& res,                                   \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    Func(res.B(), D1.B(), D2.B());                                             \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> Func                                   \
(                                                                              \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        meshDirection<ReturnType,MeshType>::New                                \
        (                                                                      \
            D1.fvMsh(),                                                        \
            D1.levelNum(),                                                     \
            D1.directionNum()                                                  \
        );                                                                     \
    Func(tRes.ref(), D1, D2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> Func                                   \
(                                                                              \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const tmp<meshDirection<Type2,MeshType>>& tD2                              \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmp<ReturnType,Type2,MeshType>::New(tD2);                \
    Func(tRes.ref(), D1, tD2());                                               \
    if (tD2.isTmp()) tD2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> Func                                   \
(                                                                              \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                             \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmp<ReturnType,Type1,MeshType>::New(tD1);                \
    Func(tRes.ref(), tD1(), D2);                                               \
    if (tD1.isTmp()) tD1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> Func                                   \
(                                                                              \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                             \
    const tmp<meshDirection<Type2,MeshType>>& tD2                              \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmpTmp<ReturnType,Type1,Type1,Type2,MeshType>::          \
        New(tD1,tD2);                                                          \
    Func(tRes.ref(), tD1(), tD2());                                            \
    if (tD1.isTmp()) tD1.clear();                                              \
    if (tD2.isTmp()) tD2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshDirection<ReturnType,MeshType>& res,                                   \
    const Type1& s1,                                                           \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    Func(res.B(), s1, D2.B());                                                 \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> Func                                   \
(                                                                              \
    const Type1& s1,                                                           \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        meshDirection<ReturnType,MeshType>::New                                \
        (                                                                      \
            D2.fvMsh(),                                                        \
            D2.levelNum(),                                                     \
            D2.directionNum()                                                  \
        );                                                                     \
    Func(tRes.ref(), s1, D2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> Func                                   \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<meshDirection<Type2,MeshType>>& tD2                              \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmp<ReturnType,Type2,MeshType>::New(tD2);                \
    Func(tRes.ref(), s1, tD2());                                               \
    if (tD2.isTmp()) tD2.clear();                                              \
    return tRes;                                                               \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    meshDirection<ReturnType,MeshType>& res,                                   \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    Func(res.B(), D1.B(), s2);                                                 \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> Func                                   \
(                                                                              \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        meshDirection<ReturnType,MeshType>::New                                \
        (                                                                      \
            D1.fvMsh(),                                                        \
            D1.levelNum(),                                                     \
            D1.directionNum()                                                  \
        );                                                                     \
    Func(tRes.ref(), D1, s2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> Func                                   \
(                                                                              \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                             \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmp<ReturnType,Type1,MeshType>::New(tD1);                \
    Func(tRes.ref(), tD1(), s2);                                               \
    if (tD1.isTmp()) tD1.clear();                                              \
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
    meshDirection<ReturnType,MeshType>& res,                                   \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    OpFunc(res.B(), D1.B(), D2.B());                                           \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> operator Op                            \
(                                                                              \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        meshDirection<ReturnType,MeshType>::New                                \
        (                                                                      \
            D1.fvMsh(),                                                        \
            D1.levelNum(),                                                     \
            D1.directionNum()                                                  \
        );                                                                     \
    OpFunc(tRes.ref(), D1, D2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> operator Op                            \
(                                                                              \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const tmp<meshDirection<Type2,MeshType>>& tD2                              \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmp<ReturnType,Type2,MeshType>::New(tD2);                \
    OpFunc(tRes.ref(), D1, tD2());                                             \
    if (tD2.isTmp()) tD2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> operator Op                            \
(                                                                              \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                             \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmp<ReturnType,Type1,MeshType>::New(tD1);                \
    OpFunc(tRes.ref(), tD1(), D2);                                             \
    if (tD1.isTmp()) tD1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> operator Op                            \
(                                                                              \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                             \
    const tmp<meshDirection<Type2,MeshType>>& tD2                              \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmpTmp<ReturnType,Type1,Type1,Type2,MeshType>::          \
        New(tD1, tD2);                                                         \
    OpFunc(tRes.ref(), tD1(), tD2());                                          \
    if (tD1.isTmp()) tD1.clear();                                              \
    if (tD2.isTmp()) tD2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshDirection<ReturnType,MeshType>& res,                                   \
    const Type1& s1,                                                           \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    OpFunc(res.B(), s1, D2.B());                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> operator Op                            \
(                                                                              \
    const Type1& s1,                                                           \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        meshDirection<ReturnType,MeshType>::New                                \
        (                                                                      \
            D2.fvMsh(),                                                        \
            D2.levelNum(),                                                     \
            D2.directionNum()                                                  \
        );                                                                     \
    OpFunc(tRes.ref(), s1, D2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> operator Op                            \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<meshDirection<Type2,MeshType>>& tD2                              \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmp<ReturnType,Type2,MeshType>::New(tD2);                \
    OpFunc(tRes.ref(), s1, tD2());                                             \
    if (tD2.isTmp()) tD2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    meshDirection<ReturnType,MeshType>& res,                                   \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    OpFunc(res.B(), D1.B(), s2);                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> operator Op                            \
(                                                                              \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        meshDirection<ReturnType,MeshType>::New                                \
        (                                                                      \
            D1.fvMsh(),                                                        \
            D1.levelNum(),                                                     \
            D1.directionNum()                                                  \
        );                                                                     \
    OpFunc(tRes.ref(), D1, s2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<meshDirection<ReturnType,MeshType>> operator Op                            \
(                                                                              \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                             \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                             \
        reuseDirectionTmp<ReturnType,Type1,MeshType>::New(tD1);                \
    OpFunc(tRes.ref(), tD1(), s2);                                             \
    if (tD1.isTmp()) tD1.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)
