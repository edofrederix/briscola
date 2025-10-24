#include "blockM.H"
#include "blockReuseFunctions.H"

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func(block<ReturnType>& res, const block<Type>& f)                        \
{                                                                              \
    BFOR_ALL_F_OP_FUNC_F(ReturnType, res, =, ::Foam::Func, Type, f)            \
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
    BFOR_ALL_F_OP_OP_F(ReturnType, res, =, Op, Type, f)                        \
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
    BFOR_ALL_F_OP_FUNC_F_F                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, f1, Type2, f2                 \
    )                                                                          \
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
    BFOR_ALL_F_OP_FUNC_S_F                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, s1, Type2, f2                 \
    )                                                                          \
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
    BFOR_ALL_F_OP_FUNC_F_S                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, f1, Type2, s2                 \
    )                                                                          \
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
    BFOR_ALL_F_OP_F_OP_F(ReturnType, res, =, Type1, f1, Op, Type2, f2)         \
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
    BFOR_ALL_F_OP_S_OP_F(ReturnType, res, =, Type1, s1, Op, Type2, f2)         \
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
    BFOR_ALL_F_OP_F_OP_S(ReturnType, res, =, Type1, f1, Op, Type2, s2)         \
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

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type1, Type2>::type>& res,                          \
    const block<Type1>& f1,                                                    \
    const block<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type ReturnType;                   \
    BFOR_ALL_F_OP_F_OP_F(ReturnType, res, =, Type1, f1, Op, Type2, f2)         \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const block<Type1>& f1, const block<Type2>& f2)                    \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType>> tRes(new block<productType>(f1.shape()));          \
    OpFunc(tRes.ref(), f1, f2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const block<Type1>& f1, const tmp<block<Type2>>& tf2)              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType>> tRes = reuseTmp<productType, Type2>::New(tf2);     \
    OpFunc(tRes.ref(), f1, tf2());                                             \
    if (tf2.isTmp()) tf2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<block<Type1>>& tf1, const block<Type2>& f2)              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType>> tRes = reuseTmp<productType, Type1>::New(tf1);     \
    OpFunc(tRes.ref(), tf1(), f2);                                             \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<block<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<block<Type1>>& tf1, const tmp<block<Type2>>& tf2)        \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<block<productType>> tRes =                                             \
        reuseTmpTmp<productType, Type1, Type1, Type2>::New(tf1, tf2);          \
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
    direction nCmpt                                                            \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type, Form>::type>& res,                            \
    const block<Type>& f1,                                                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type ReturnType;                     \
    checkBlocks(res, f1, "res = f1 " #Op " s");                                \
    List_ACCESS(ReturnType, res, resP);                                        \
    List_CONST_ACCESS(Type, f1, f1P);                                          \
    List_FOR_ALL(res, i)                                                       \
        List_ELEM(res, resP, i) =                                              \
            List_ELEM(f1, f1P, i) Op static_cast<const Form&>(vs);             \
    List_END_FOR_ALL                                                           \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op(const block<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs)     \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType>> tRes(new block<productType>(f1.shape()));          \
    OpFunc(tRes.ref(), f1, static_cast<const Form&>(vs));                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const tmp<block<Type>>& tf1,                                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType>> tRes = reuseTmp<productType, Type>::New(tf1);      \
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
    direction nCmpt                                                            \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Form, Type>::type>& res,                            \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const block<Type>& f1                                                      \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type ReturnType;                     \
    checkBlocks(res, f1, "res = s " #Op " f1");                                \
    List_ACCESS(ReturnType, res, resP);                                        \
    List_CONST_ACCESS(Type, f1, f1P);                                          \
    List_FOR_ALL(res, i)                                                       \
        List_ELEM(res, resP, i) =                                              \
            static_cast<const Form&>(vs) Op List_ELEM(f1, f1P, i);             \
    List_END_FOR_ALL                                                           \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const block<Type>& f1)     \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType>> tRes(new block<productType>(f1.shape()));          \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), f1);                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<block<Type>>& tf1                                                \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType>> tRes = reuseTmp<productType, Type>::New(tf1);      \
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
    direction nCmpt                                                            \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Type, Form>::type>& res,                            \
    const block<Type>& f1,                                                     \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type ReturnType;                     \
    checkBlocks(res, f1, "res = f1 " #Op " s");                                \
    List_ACCESS(ReturnType, res, resP);                                        \
    List_CONST_ACCESS(Type, f1, f1P);                                          \
    List_FOR_ALL(res, i)                                                       \
        List_ELEM(res, resP, i) =                                              \
            List_ELEM(f1, f1P, i) Op static_cast<const Form&>(vs);             \
    List_END_FOR_ALL                                                           \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op(const block<Type>& f1, const CellSpace<Form,Cmpt,nCmpt>& vs)       \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType>> tRes(new block<productType>(f1.shape()));          \
    OpFunc(tRes.ref(), f1, static_cast<const Form&>(vs));                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Type, Form>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const tmp<block<Type>>& tf1,                                               \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<block<productType>> tRes = reuseTmp<productType, Type>::New(tf1);      \
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
    direction nCmpt                                                            \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    block<typename product<Form, Type>::type>& res,                            \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const block<Type>& f1                                                      \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type ReturnType;                     \
    checkBlocks(res, f1, "res = s " #Op " f1");                                \
    List_ACCESS(ReturnType, res, resP);                                        \
    List_CONST_ACCESS(Type, f1, f1P);                                          \
    List_FOR_ALL(res, i)                                                       \
        List_ELEM(res, resP, i) =                                              \
            static_cast<const Form&>(vs) Op List_ELEM(f1, f1P, i);             \
    List_END_FOR_ALL                                                           \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op(const CellSpace<Form,Cmpt,nCmpt>& vs, const block<Type>& f1)       \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType>> tRes(new block<productType>(f1.shape()));          \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), f1);                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt                                                            \
>                                                                              \
tmp<block<typename product<Form, Type>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const tmp<block<Type>>& tf1                                                \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<block<productType>> tRes = reuseTmp<productType, Type>::New(tf1);      \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), tf1());                   \
    if (tf1.isTmp()) tf1.clear();                                              \
    return tRes;                                                               \
}
