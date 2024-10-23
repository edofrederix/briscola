#include "meshDirectionReuseFunctions.H"

#define UNARY_FUNCTION(ReturnType, Type, Func)                                  \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshDirection<ReturnType,MeshType>& res,                                    \
    const meshDirection<Type,MeshType>& D                                       \
)                                                                               \
{                                                                               \
    Func(res.B(), D.B());                                                       \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>>                                         \
Func(const meshDirection<Type,MeshType>& D)                                     \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes                                \
    (                                                                           \
        new meshDirection<ReturnType,MeshType>                                  \
        (                                                                       \
            D.fvMsh(),                                                          \
            D.levelNum(),                                                       \
            D.directionNum()                                                    \
        )                                                                       \
    );                                                                          \
    Func(tRes.ref(), D);                                                        \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>>                                         \
Func(const tmp<meshDirection<Type,MeshType>>& tD)                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmp<ReturnType,Type,MeshType>::New(tD);                         \
    Func(tRes.ref(), tD());                                                     \
    if (tD.isTmp()) tD.clear();                                                 \
    return tRes;                                                                \
}

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                            \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshDirection<ReturnType,MeshType>& res,                                    \
    const meshDirection<Type,MeshType>& D                                       \
)                                                                               \
{                                                                               \
    OpFunc(res.B(), D.B());                                                     \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>>                                         \
operator Op(const meshDirection<Type,MeshType>& D)                              \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes                                \
    (                                                                           \
        new meshDirection<ReturnType,MeshType>                                  \
        (                                                                       \
            D.fvMsh(),                                                          \
            D.levelNum(),                                                       \
            D.directionNum()                                                    \
        )                                                                       \
    );                                                                          \
    OpFunc(tRes.ref(), D);                                                      \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>>                                         \
operator Op(const tmp<meshDirection<Type,MeshType>>& tD)                        \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmp<ReturnType,Type,MeshType>::New(tD);                         \
    OpFunc(tRes.ref(), tD());                                                   \
    if (tD.isTmp()) tD.clear();                                                 \
    return tRes;                                                                \
}

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                         \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshDirection<ReturnType,MeshType>& res,                                    \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    Func(res.B(), D1.B(), D2.B());                                              \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes                                \
    (                                                                           \
        new meshDirection<ReturnType,MeshType>                                  \
        (                                                                       \
            D1.fvMsh(),                                                         \
            D1.levelNum(),                                                      \
            D1.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    Func(tRes.ref(), D1, D2);                                                   \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const tmp<meshDirection<Type2,MeshType>>& tD2                               \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmp<ReturnType,Type2,MeshType>::New(tD2);                       \
    Func(tRes.ref(), D1, tD2());                                                \
    if (tD2.isTmp()) tD2.clear();                                               \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                              \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmp<ReturnType,Type1,MeshType>::New(tD1);                       \
    Func(tRes.ref(), tD1(), D2);                                                \
    if (tD1.isTmp()) tD1.clear();                                               \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                              \
    const tmp<meshDirection<Type2,MeshType>>& tD2                               \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmpTmp<ReturnType,Type1,Type1,Type2,MeshType>::New(tD1,tD2);    \
    Func(tRes.ref(), tD1(), tD2());                                             \
    if (tD1.isTmp()) tD1.clear();                                               \
    if (tD2.isTmp()) tD2.clear();                                               \
    return tRes;                                                                \
}

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshDirection<ReturnType,MeshType>& res,                                    \
    const Type1& s1,                                                            \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    Func(res.B(), s1, D2.B());                                                  \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const Type1& s1,                                                            \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes                                \
    (                                                                           \
        new meshDirection<ReturnType,MeshType>                                  \
        (                                                                       \
            D2.fvMsh(),                                                         \
            D2.levelNum(),                                                      \
            D2.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    Func(tRes.ref(), s1, D2);                                                   \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const Type1& s1,                                                            \
    const tmp<meshDirection<Type2,MeshType>>& tD2                               \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmp<ReturnType,Type2,MeshType>::New(tD2);                       \
    Func(tRes.ref(), s1, tD2());                                                \
    if (tD2.isTmp()) tD2.clear();                                               \
    return tRes;                                                                \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshDirection<ReturnType,MeshType>& res,                                    \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    Func(res.B(), D1.B(), s2);                                                  \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes                                \
    (                                                                           \
        new meshDirection<ReturnType,MeshType>                                  \
        (                                                                       \
            D1.fvMsh(),                                                         \
            D1.levelNum(),                                                      \
            D1.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    Func(tRes.ref(), D1, s2);                                                   \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                              \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmp<ReturnType,Type1,MeshType>::New(tD1);                       \
    Func(tRes.ref(), tD1(), s2);                                                \
    if (tD1.isTmp()) tD1.clear();                                               \
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
    meshDirection<ReturnType,MeshType>& res,                                    \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    OpFunc(res.B(), D1.B(), D2.B());                                            \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes                                \
    (                                                                           \
        new meshDirection<ReturnType,MeshType>                                  \
        (                                                                       \
            D1.fvMsh(),                                                         \
            D1.levelNum(),                                                      \
            D1.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    OpFunc(tRes.ref(), D1, D2);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const tmp<meshDirection<Type2,MeshType>>& tD2                               \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmp<ReturnType,Type2,MeshType>::New(tD2);                       \
    OpFunc(tRes.ref(), D1, tD2());                                              \
    if (tD2.isTmp()) tD2.clear();                                               \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                              \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmp<ReturnType,Type1,MeshType>::New(tD1);                       \
    OpFunc(tRes.ref(), tD1(), D2);                                              \
    if (tD1.isTmp()) tD1.clear();                                               \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                              \
    const tmp<meshDirection<Type2,MeshType>>& tD2                               \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmpTmp<ReturnType,Type1,Type1,Type2,MeshType>::New(tD1, tD2);   \
    OpFunc(tRes.ref(), tD1(), tD2());                                           \
    if (tD1.isTmp()) tD1.clear();                                               \
    if (tD2.isTmp()) tD2.clear();                                               \
    return tRes;                                                                \
}

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshDirection<ReturnType,MeshType>& res,                                    \
    const Type1& s1,                                                            \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    OpFunc(res.B(), s1, D2.B());                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const Type1& s1,                                                            \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes                                \
    (                                                                           \
        new meshDirection<ReturnType,MeshType>                                  \
        (                                                                       \
            D2.fvMsh(),                                                         \
            D2.levelNum(),                                                      \
            D2.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    OpFunc(tRes.ref(), s1, D2);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const Type1& s1,                                                            \
    const tmp<meshDirection<Type2,MeshType>>& tD2                               \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmp<ReturnType,Type2,MeshType>::New(tD2);                       \
    OpFunc(tRes.ref(), s1, tD2());                                              \
    if (tD2.isTmp()) tD2.clear();                                               \
    return tRes;                                                                \
}

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshDirection<ReturnType,MeshType>& res,                                    \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    OpFunc(res.B(), D1.B(), s2);                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes                                \
    (                                                                           \
        new meshDirection<ReturnType,MeshType>                                  \
        (                                                                       \
            D1.fvMsh(),                                                         \
            D1.levelNum(),                                                      \
            D1.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    OpFunc(tRes.ref(), D1, s2);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                              \
    const Type2& s2                                                             \
)                                                                               \
{                                                                               \
    tmp<meshDirection<ReturnType,MeshType>> tRes =                              \
        reuseDirTmp<ReturnType,Type1,MeshType>::New(tD1);                       \
    OpFunc(tRes.ref(), tD1(), s2);                                              \
    if (tD1.isTmp()) tD1.clear();                                               \
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
    meshDirection<typename product<Type1, Type2>::type,MeshType>& res,          \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    OpFunc(res.B(), D1.B(), D2.B());                                            \
}                                                                               \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>               \
operator Op                                                                     \
(                                                                               \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    typedef typename product<Type1, Type2>::type productType;                   \
    tmp<meshDirection<productType,MeshType>> tRes                               \
    (                                                                           \
        new meshDirection<productType,MeshType>                                 \
        (                                                                       \
            D1.fvMsh(),                                                         \
            D1.levelNum(),                                                      \
            D1.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    OpFunc(tRes.ref(), D1, D2);                                                 \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>               \
operator Op                                                                     \
(                                                                               \
    const meshDirection<Type1,MeshType>& D1,                                    \
    const tmp<meshDirection<Type2,MeshType>>& tD2                               \
)                                                                               \
{                                                                               \
    typedef typename product<Type1, Type2>::type productType;                   \
    tmp<meshDirection<productType,MeshType>> tRes =                             \
        reuseDirTmp<productType,Type2,MeshType>::New(tD2);                      \
    OpFunc(tRes.ref(), D1, tD2());                                              \
    if (tD2.isTmp()) tD2.clear();                                               \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>               \
operator Op                                                                     \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                              \
    const meshDirection<Type2,MeshType>& D2                                     \
)                                                                               \
{                                                                               \
    typedef typename product<Type1, Type2>::type productType;                   \
    tmp<meshDirection<productType,MeshType>> tRes =                             \
        reuseDirTmp<productType,Type1,MeshType>::New(tD1);                      \
    OpFunc(tRes.ref(), tD1(), D2);                                              \
    if (tD1.isTmp()) tD1.clear();                                               \
    return tRes;                                                                \
}                                                                               \
                                                                                \
template<class Type1, class Type2, class MeshType>                              \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>               \
operator Op                                                                     \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                              \
    const tmp<meshDirection<Type2,MeshType>>& tD2                               \
)                                                                               \
{                                                                               \
    typedef typename product<Type1, Type2>::type productType;                   \
    tmp<meshDirection<productType,MeshType>> tRes =                             \
        reuseDirTmpTmp<productType,Type1,Type1,Type2,MeshType>::New(tD1, tD2);  \
    OpFunc(tRes.ref(), tD1(), tD2());                                           \
    if (tD1.isTmp()) tD1.clear();                                               \
    if (tD2.isTmp()) tD2.clear();                                               \
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
    meshDirection<typename product<Type, Form>::type,MeshType>& res,            \
    const meshDirection<Type,MeshType>& D1,                                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                      \
)                                                                               \
{                                                                               \
    OpFunc(res.B(), D1.B(), vs);                                                \
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
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                 \
operator Op                                                                     \
(                                                                               \
    const meshDirection<Type,MeshType>& D1,                                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                      \
)                                                                               \
{                                                                               \
    typedef typename product<Type, Form>::type productType;                     \
    tmp<meshDirection<productType,MeshType>> tRes                               \
    (                                                                           \
        new meshDirection<productType,MeshType>                                 \
        (                                                                       \
            D1.fvMsh(),                                                         \
            D1.levelNum(),                                                      \
            D1.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    OpFunc(tRes.ref(), D1, static_cast<const Form&>(vs));                       \
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
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                 \
operator Op                                                                     \
(                                                                               \
    const tmp<meshDirection<Type,MeshType>>& tD1,                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                      \
)                                                                               \
{                                                                               \
    typedef typename product<Type, Form>::type productType;                     \
    tmp<meshDirection<productType,MeshType>> tRes =                             \
        reuseDirTmp<productType,Type,MeshType>::New(tD1);                       \
    OpFunc(tRes.ref(), tD1(), static_cast<const Form&>(vs));                    \
    if (tD1.isTmp()) tD1.clear();                                               \
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
    meshDirection<typename product<Form, Type>::type,MeshType>& res,            \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                     \
    const meshDirection<Type,MeshType>& D1                                      \
)                                                                               \
{                                                                               \
    OpFunc(res.B(), vs, D1.B());                                                \
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
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                 \
operator Op                                                                     \
(                                                                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                     \
    const meshDirection<Type,MeshType>& D1                                      \
)                                                                               \
{                                                                               \
    typedef typename product<Form, Type>::type productType;                     \
    tmp<meshDirection<productType,MeshType>> tRes                               \
    (                                                                           \
        new meshDirection<productType,MeshType>                                 \
        (                                                                       \
            D1.fvMsh(),                                                         \
            D1.levelNum(),                                                      \
            D1.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), D1);                       \
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
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                 \
operator Op                                                                     \
(                                                                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                     \
    const tmp<meshDirection<Type,MeshType>>& tD1                                \
)                                                                               \
{                                                                               \
    typedef typename product<Form, Type>::type productType;                     \
    tmp<meshDirection<productType,MeshType>> tRes =                             \
        reuseDirTmp<productType,Type,MeshType>::New(tD1);                       \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), tD1());                    \
    if (tD1.isTmp()) tD1.clear();                                               \
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
    meshDirection<typename product<Type, Form>::type,MeshType>& res,            \
    const meshDirection<Type,MeshType>& D1,                                     \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                        \
)                                                                               \
{                                                                               \
    OpFunc(res.B(), D1.B(), vs);                                                \
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
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                 \
operator Op                                                                     \
(                                                                               \
    const meshDirection<Type,MeshType>& D1,                                     \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                        \
)                                                                               \
{                                                                               \
    typedef typename product<Type, Form>::type productType;                     \
    tmp<meshDirection<productType,MeshType>> tRes                               \
    (                                                                           \
        new meshDirection<productType,MeshType>                                 \
        (                                                                       \
            D1.fvMsh(),                                                         \
            D1.levelNum(),                                                      \
            D1.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    OpFunc(tRes.ref(), D1, static_cast<const Form&>(vs));                       \
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
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                 \
operator Op                                                                     \
(                                                                               \
    const tmp<meshDirection<Type,MeshType>>& tD1,                               \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                        \
)                                                                               \
{                                                                               \
    typedef typename product<Type, Form>::type productType;                     \
    tmp<meshDirection<productType,MeshType>> tRes =                             \
        reuseDirTmp<productType,Type,MeshType>::New(tD1);                       \
    OpFunc(tRes.ref(), tD1(), static_cast<const Form&>(vs));                    \
    if (tD1.isTmp()) tD1.clear();                                               \
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
    meshDirection<typename product<Form, Type>::type,MeshType>& res,            \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                       \
    const meshDirection<Type,MeshType>& D1                                      \
)                                                                               \
{                                                                               \
    OpFunc(res.B(), vs, D1.B());                                                \
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
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                 \
operator Op                                                                     \
(                                                                               \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                       \
    const meshDirection<Type,MeshType>& D1                                      \
)                                                                               \
{                                                                               \
    typedef typename product<Form, Type>::type productType;                     \
    tmp<meshDirection<productType,MeshType>> tRes                               \
    (                                                                           \
        new meshDirection<productType,MeshType>                                 \
        (                                                                       \
            D1.fvMsh(),                                                         \
            D1.levelNum(),                                                      \
            D1.directionNum()                                                   \
        )                                                                       \
    );                                                                          \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), D1);                       \
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
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                 \
operator Op                                                                     \
(                                                                               \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                       \
    const tmp<meshDirection<Type,MeshType>>& tD1                                \
)                                                                               \
{                                                                               \
    typedef typename product<Form, Type>::type productType;                     \
    tmp<meshDirection<productType,MeshType>> tRes =                             \
        reuseDirTmp<productType,Type,MeshType>::New(tD1);                       \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), tD1());                    \
    if (tD1.isTmp()) tD1.clear();                                               \
    return tRes;                                                                \
}
