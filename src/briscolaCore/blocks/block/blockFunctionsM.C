#include "blockM.H"
#include "blockReuseFunctions.H"

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func(block<ReturnType,P>& res, const block<Type,P>& f)                    \
{                                                                              \
    checkBlocks<ReturnType,Type,P>(res,f,#Func);                               \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = ::Foam::Func(f(i));                                           \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func(const block<Type,P>& f)                          \
{                                                                              \
    tmp<block<ReturnType,P>> tRes(new block<ReturnType,P>(f.shape()));         \
    Func(tRes.ref(), f);                                                       \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func(const tmp<block<Type,P>>& tf)                    \
{                                                                              \
    tmp<block<ReturnType,P>> tRes = reuseTmp<ReturnType, Type, P>::New(tf);    \
    Func(tRes.ref(), tf());                                                    \
    if (tf.isTmp()) tf.clear();                                                \
    return tRes;                                                               \
}

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                           \
                                                                               \
TEMPLATE                                                                       \
void OpFunc(block<ReturnType,P>& res, const block<Type,P>& f)                  \
{                                                                              \
    checkBlocks<ReturnType,Type,P>(res,f,#OpFunc);                             \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = Op f(i);                                                      \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op(const block<Type,P>& f)                   \
{                                                                              \
    tmp<block<ReturnType,P>> tRes(new block<ReturnType,P>(f.shape()));         \
    OpFunc(tRes.ref(), f);                                                     \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op(const tmp<block<Type,P>>& tf)             \
{                                                                              \
    tmp<block<ReturnType,P>> tRes = reuseTmp<ReturnType, Type, P>::New(tf);    \
    OpFunc(tRes.ref(), tf());                                                  \
    if (tf.isTmp()) tf.clear();                                                \
    return tRes;                                                               \
}

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    block<ReturnType,P>& res,                                                  \
    const block<Type1,P>& f1,                                                  \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    checkBlocks<ReturnType,Type1,Type2,P>(res,f1,f2,#Func);                    \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = ::Foam::Func(f1(i),f2(i));                                    \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const block<Type1,P>& f1,                                                  \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes(new block<ReturnType,P>(f1.shape()));        \
    Func(tRes.ref(), f1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const block<Type1,P>& f1,                                                  \
    const tmp<block<Type2,P>>& tf2                                             \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes = reuseTmp<ReturnType, Type2, P>::New(tf2);  \
    Func(tRes.ref(), f1, tf2());                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const tmp<block<Type1,P>>& tf1,                                            \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes = reuseTmp<ReturnType, Type1, P>::New(tf1);  \
    Func(tRes.ref(), tf1(), f2);                                               \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const tmp<block<Type1,P>>& tf1,                                            \
    const tmp<block<Type2,P>>& tf2                                             \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes =                                            \
        reuseTmpTmp<ReturnType, Type1, Type1, Type2, P>::New(tf1, tf2);        \
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
    block<ReturnType,P>& res,                                                  \
    const Type1& s1,                                                           \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    checkBlocks<ReturnType,Type2,P>(res,f2,#Func);                             \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = ::Foam::Func(s1,f2(i));                                       \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const Type1& s1,                                                           \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes(new block<ReturnType,P>(f2.shape()));        \
    Func(tRes.ref(), s1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<block<Type2,P>>& tf2                                             \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes = reuseTmp<ReturnType, Type2, P>::New(tf2);  \
    Func(tRes.ref(), s1, tf2());                                               \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    block<ReturnType,P>& res,                                                  \
    const block<Type1,P>& f1,                                                  \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    checkBlocks<ReturnType,Type1,P>(res,f1,#Func);                             \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = ::Foam::Func(f1(i),s2);                                       \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const block<Type1,P>& f1,                                                  \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes(new block<ReturnType,P>(f1.shape()));        \
    Func(tRes.ref(), f1, s2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const tmp<block<Type1,P>>& tf1,                                            \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes = reuseTmp<ReturnType, Type1, P>::New(tf1);  \
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
    block<ReturnType,P>& res,                                                  \
    const block<Type1,P>& f1,                                                  \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    checkBlocks<ReturnType,Type1,Type2,P>(res,f1,f2,#OpFunc);                  \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = f1(i) Op f2(i);                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const block<Type1,P>& f1,                                                  \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes(new block<ReturnType,P>(f1.shape()));        \
    OpFunc(tRes.ref(), f1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const block<Type1,P>& f1,                                                  \
    const tmp<block<Type2,P>>& tf2                                             \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes = reuseTmp<ReturnType, Type2, P>::New(tf2);  \
    OpFunc(tRes.ref(), f1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const tmp<block<Type1,P>>& tf1,                                            \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes = reuseTmp<ReturnType, Type1, P>::New(tf1);  \
    OpFunc(tRes.ref(), tf1(), f2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const tmp<block<Type1,P>>& tf1,                                            \
    const tmp<block<Type2,P>>& tf2                                             \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes =                                            \
        reuseTmpTmp<ReturnType, Type1, Type1, Type2, P>::New(tf1, tf2);        \
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
    block<ReturnType,P>& res,                                                  \
    const Type1& s1,                                                           \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    checkBlocks<ReturnType,Type2,P>(res,f2,#OpFunc);                           \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = s1 Op f2(i);                                                  \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const Type1& s1,                                                           \
    const block<Type2,P>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes(new block<ReturnType,P>(f2.shape()));        \
    OpFunc(tRes.ref(), s1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<block<Type2,P>>& tf2                                             \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes = reuseTmp<ReturnType, Type2, P>::New(tf2);  \
    OpFunc(tRes.ref(), s1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    block<ReturnType,P>& res,                                                  \
    const block<Type1,P>& f1,                                                  \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    checkBlocks<ReturnType,Type1,P>(res,f1,#OpFunc);                           \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
        res(i) = f1(i) Op s2;                                                  \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const block<Type1,P>& f1,                                                  \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes(new block<ReturnType,P>(f1.shape()));        \
    OpFunc(tRes.ref(), f1, s2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const tmp<block<Type1,P>>& tf1,                                            \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<block<ReturnType,P>> tRes = reuseTmp<ReturnType, Type1, P>::New(tf1);  \
    OpFunc(tRes.ref(), tf1(), s2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)
