#include "blockM.H"
#include "blockReuseFunctions.H"

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func(block<ReturnType>& res, const block<Type>& f)                        \
{                                                                              \
    checkMatrices(res,f,#Func);                                                \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = ::Foam::Func(f(i));                                           \
    }                                                                          \
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
    tf.clear();                                                                \
    return tRes;                                                               \
}

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                           \
                                                                               \
TEMPLATE                                                                       \
void OpFunc(block<ReturnType>& res, const block<Type>& f)                      \
{                                                                              \
    checkMatrices(res,f,#OpFunc);                                              \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = Op f(i);                                                      \
    }                                                                          \
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
    tf.clear();                                                                \
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
    checkMatrices(res,f1,f2,#Func);                                            \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = ::Foam::Func(f1(i),f2(i));                                    \
    }                                                                          \
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
    tf2.clear();                                                               \
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
    tf1.clear();                                                               \
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
    tf1.clear();                                                               \
    tf2.clear();                                                               \
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
    checkMatrices(res,f2,#Func);                                               \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = ::Foam::Func(s1,f2(i));                                       \
    }                                                                          \
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
    tf2.clear();                                                               \
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
    checkMatrices(res,f1,#Func);                                               \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = ::Foam::Func(f1(i),s2);                                       \
    }                                                                          \
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
    tf1.clear();                                                               \
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
    checkMatrices(res,f1,f2,#OpFunc);                                          \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = f1(i) Op f2(i);                                               \
    }                                                                          \
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
    tf2.clear();                                                               \
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
    tf1.clear();                                                               \
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
    tf1.clear();                                                               \
    tf2.clear();                                                               \
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
    checkMatrices(res,f2,#OpFunc);                                             \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = s1 Op f2(i);                                                  \
    }                                                                          \
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
    tf2.clear();                                                               \
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
    checkMatrices(res,f1,#OpFunc);                                             \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = f1(i) Op s2;                                                  \
    }                                                                          \
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
    tf1.clear();                                                               \
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
    checkMatrices(res,f1,f2,#OpFunc);                                          \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = f1(i) Op f2(i);                                               \
    }                                                                          \
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
    tf2.clear();                                                               \
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
    tf1.clear();                                                               \
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
    tf1.clear();                                                               \
    tf2.clear();                                                               \
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
    checkMatrices(res,f1,#OpFunc);                                             \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = f1(i) Op static_cast<const Form&>(vs);                        \
    }                                                                          \
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
    tf1.clear();                                                               \
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
    checkMatrices(res,f1,#OpFunc);                                             \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = static_cast<const Form&>(vs) Op f1(i);                        \
    }                                                                          \
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
    tf1.clear();                                                               \
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
    checkMatrices(res,f1,#OpFunc);                                             \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = f1(i) Op static_cast<const Form&>(vs);                        \
    }                                                                          \
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
    tf1.clear();                                                               \
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
    checkMatrices(res,f1,#OpFunc);                                             \
                                                                               \
    forAllBlockLinear(res, i)                                                  \
    {                                                                          \
        res(i) = static_cast<const Form&>(vs) Op f1(i);                        \
    }                                                                          \
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
    tf1.clear();                                                               \
    return tRes;                                                               \
}
