#define forAllBlock(BLOCK, I, J, K)                                            \
    for (int I = 0; I < (BLOCK).l(); I++)                                      \
    for (int J = 0; J < (BLOCK).m(); J++)                                      \
    for (int K = 0; K < (BLOCK).n(); K++)

#define forAllBlockLinear(BLOCK, L)                                            \
    for (int L = 0; L < (BLOCK).size(); L++)

#define UNARY_FUNCTION(ReturnType, Type1, Func)                                \
                                                                               \
TEMPLATE                                                                       \
void Func(block<ReturnType,P>&, const block<Type1,P>&);                        \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func(const block<Type1,P>&);                          \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func(const tmp<block<Type1,P>>&);

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc(block<ReturnType,P>&, const block<Type1,P>&);                      \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op(const block<Type1,P>&);                   \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op(const tmp<block<Type1,P>>&);

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    block<ReturnType,P>&,                                                      \
    const block<Type1,P>&,                                                     \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const block<Type1,P>&,                                                     \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const block<Type1,P>&,                                                     \
    const tmp<block<Type2,P>>&                                                 \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const tmp<block<Type1,P>>&,                                                \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const tmp<block<Type1,P>>&,                                                \
    const tmp<block<Type2,P>>&                                                 \
);

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    block<ReturnType,P>&,                                                      \
    const Type1&,                                                              \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const Type1&,                                                              \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const Type1&,                                                              \
    const tmp<block<Type2,P>>&                                                 \
);

#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    block<ReturnType,P>&,                                                      \
    const block<Type1,P>&,                                                     \
    const Type2&                                                               \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const block<Type1,P>&,                                                     \
    const Type2&                                                               \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> Func                                                  \
(                                                                              \
    const tmp<block<Type1,P>>&,                                                \
    const Type2&                                                               \
);

#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    block<ReturnType,P>&,                                                      \
    const block<Type1,P>&,                                                     \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const block<Type1,P>&,                                                     \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const block<Type1,P>&,                                                     \
    const tmp<block<Type2,P>>&                                                 \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const tmp<block<Type1,P>>&,                                                \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const tmp<block<Type1,P>>&,                                                \
    const tmp<block<Type2,P>>&                                                 \
);

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    block<ReturnType,P>&,                                                      \
    const Type1&,                                                              \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const Type1&,                                                              \
    const block<Type2,P>&                                                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const Type1&,                                                              \
    const tmp<block<Type2,P>>&                                                 \
);

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    block<ReturnType,P>&,                                                      \
    const block<Type1,P>&,                                                     \
    const Type2&                                                               \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const block<Type1,P>&,                                                     \
    const Type2&                                                               \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<block<ReturnType,P>> operator Op                                           \
(                                                                              \
    const tmp<block<Type1,P>>&,                                                \
    const Type2&                                                               \
);

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)
