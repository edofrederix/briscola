#include "blockM.H"
#include "blockReuseFunctions.H"

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func(block<ReturnType>& res, const block<Type>& f)                        \
{                                                                              \
    checkBlocks(res,f,#Func);                                                  \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = ::Foam::Func(f(i));                                           \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> Func(const block<Type>& f)                              \
{                                                                              \
    tmp<block<ReturnType>> tRes(new block<ReturnType>(f.shape()));             \
    Func(tRes.ref(), f);                                                       \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> Func(const tmp<block<Type>>& tf)                        \
{                                                                              \
    tmp<block<ReturnType>> tRes = reuseTmp<ReturnType, Type>::New(tf);         \
    Func(tRes.ref(), tf());                                                    \
    if (tf.isTmp()) tf.clear();                                                \
    return tRes;                                                               \
}

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                           \
                                                                               \
TEMPLATE                                                                       \
void OpFunc(block<ReturnType>& res, const block<Type>& f)                      \
{                                                                              \
    checkBlocks(res,f,#OpFunc);                                                \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = Op f(i);                                                      \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> operator Op(const block<Type>& f)                       \
{                                                                              \
    tmp<block<ReturnType>> tRes(new block<ReturnType>(f.shape()));             \
    OpFunc(tRes.ref(), f);                                                     \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> operator Op(const tmp<block<Type>>& tf)                 \
{                                                                              \
    tmp<block<ReturnType>> tRes = reuseTmp<ReturnType, Type>::New(tf);         \
    OpFunc(tRes.ref(), tf());                                                  \
    if (tf.isTmp()) tf.clear();                                                \
    return tRes;                                                               \
}

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    block<ReturnType>& res,                                                    \
    const block<Type1>& f1,                                                    \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    checkBlocks(res,f1,f2,#Func);                                              \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = ::Foam::Func(f1(i),f2(i));                                    \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> Func                                                    \
(                                                                              \
    const block<Type1>& f1,                                                    \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes(new block<ReturnType>(f1.shape()));            \
    Func(tRes.ref(), f1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> Func                                                    \
(                                                                              \
    const block<Type1>& f1,                                                    \
    const tmp<block<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes = reuseTmp<ReturnType, Type2>::New(tf2);       \
    Func(tRes.ref(), f1, tf2());                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> Func                                                    \
(                                                                              \
    const tmp<block<Type1>>& tf1,                                              \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes = reuseTmp<ReturnType, Type1>::New(tf1);       \
    Func(tRes.ref(), tf1(), f2);                                               \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> Func                                                    \
(                                                                              \
    const tmp<block<Type1>>& tf1,                                              \
    const tmp<block<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes =                                              \
        reuseTmpTmp<ReturnType, Type1, Type1, Type2>::New(tf1, tf2);           \
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
    block<ReturnType>& res,                                                    \
    const Type1& s1,                                                           \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    checkBlocks(res,f2,#Func);                                                 \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = ::Foam::Func(s1,f2(i));                                       \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> Func                                                    \
(                                                                              \
    const Type1& s1,                                                           \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes(new block<ReturnType>(f2.shape()));            \
    Func(tRes.ref(), s1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> Func                                                    \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<block<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes = reuseTmp<ReturnType, Type2>::New(tf2);       \
    Func(tRes.ref(), s1, tf2());                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    block<ReturnType>& res,                                                    \
    const block<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    checkBlocks(res,f1,#Func);                                                 \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = ::Foam::Func(f1(i),s2);                                       \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> Func                                                    \
(                                                                              \
    const block<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes(new block<ReturnType>(f1.shape()));            \
    Func(tRes.ref(), f1, s2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> Func                                                    \
(                                                                              \
    const tmp<block<Type1>>& tf1,                                              \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes = reuseTmp<ReturnType, Type1>::New(tf1);       \
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
    block<ReturnType>& res,                                                    \
    const block<Type1>& f1,                                                    \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    checkBlocks(res,f1,f2,#OpFunc);                                            \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = f1(i) Op f2(i);                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> operator Op                                             \
(                                                                              \
    const block<Type1>& f1,                                                    \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes(new block<ReturnType>(f1.shape()));            \
    OpFunc(tRes.ref(), f1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> operator Op                                             \
(                                                                              \
    const block<Type1>& f1,                                                    \
    const tmp<block<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes = reuseTmp<ReturnType, Type2>::New(tf2);       \
    OpFunc(tRes.ref(), f1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> operator Op                                             \
(                                                                              \
    const tmp<block<Type1>>& tf1,                                              \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes = reuseTmp<ReturnType, Type1>::New(tf1);       \
    OpFunc(tRes.ref(), tf1(), f2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> operator Op                                             \
(                                                                              \
    const tmp<block<Type1>>& tf1,                                              \
    const tmp<block<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes =                                              \
        reuseTmpTmp<ReturnType, Type1, Type1, Type2>::New(tf1, tf2);           \
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
    block<ReturnType>& res,                                                    \
    const Type1& s1,                                                           \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    checkBlocks(res,f2,#OpFunc);                                               \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = s1 Op f2(i);                                                  \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> operator Op                                             \
(                                                                              \
    const Type1& s1,                                                           \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes(new block<ReturnType>(f2.shape()));            \
    OpFunc(tRes.ref(), s1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> operator Op                                             \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<block<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes = reuseTmp<ReturnType, Type2>::New(tf2);       \
    OpFunc(tRes.ref(), s1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    block<ReturnType>& res,                                                    \
    const block<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    checkBlocks(res,f1,#OpFunc);                                               \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = f1(i) Op s2;                                                  \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> operator Op                                             \
(                                                                              \
    const block<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes(new block<ReturnType>(f1.shape()));            \
    OpFunc(tRes.ref(), f1, s2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType>> operator Op                                             \
(                                                                              \
    const tmp<block<Type1>>& tf1,                                              \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType>> tRes = reuseTmp<ReturnType, Type1>::New(tf1);       \
    OpFunc(tRes.ref(), tf1(), s2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)
