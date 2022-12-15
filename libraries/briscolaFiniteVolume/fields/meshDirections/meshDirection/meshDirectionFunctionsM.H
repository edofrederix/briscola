#define UNARY_FUNCTION(ReturnType, Type1, Func)                                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshDirection<ReturnType,MeshType>&,                                        \
    const meshDirection<Type1,MeshType>&                                        \
);                                                                              \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>>                                         \
Func(const meshDirection<Type1,MeshType>&);                                     \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>>                                         \
Func(const tmp<meshDirection<Type1,MeshType>>&);

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshDirection<ReturnType,MeshType>&,                                        \
    const meshDirection<Type1,MeshType>&                                        \
);                                                                              \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>>                                         \
operator Op(const meshDirection<Type1,MeshType>&);                              \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>>                                         \
operator Op(const tmp<meshDirection<Type1,MeshType>>&);

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                         \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshDirection<ReturnType,MeshType>&,                                        \
    const meshDirection<Type1,MeshType>&,                                       \
    const meshDirection<Type2,MeshType>&                                        \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const meshDirection<Type1,MeshType>&,                                       \
    const meshDirection<Type2,MeshType>&                                        \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const meshDirection<Type1,MeshType>&,                                       \
    const tmp<meshDirection<Type2,MeshType>>&                                   \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>&,                                  \
    const meshDirection<Type2,MeshType>&                                        \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>&,                                  \
    const tmp<meshDirection<Type2,MeshType>>&                                   \
);

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshDirection<ReturnType,MeshType>&,                                        \
    const Type1&,                                                               \
    const meshDirection<Type2,MeshType>&                                        \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const Type1&,                                                               \
    const meshDirection<Type2,MeshType>&                                        \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const Type1&,                                                               \
    const tmp<meshDirection<Type2,MeshType>>&                                   \
);

#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshDirection<ReturnType,MeshType>&,                                        \
    const meshDirection<Type1,MeshType>&,                                       \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const meshDirection<Type1,MeshType>&,                                       \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> Func                                    \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>&,                                  \
    const Type2&                                                                \
);

#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                     \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                   \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshDirection<ReturnType,MeshType>&,                                        \
    const meshDirection<Type1,MeshType>&,                                       \
    const meshDirection<Type2,MeshType>&                                        \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const meshDirection<Type1,MeshType>&,                                       \
    const meshDirection<Type2,MeshType>&                                        \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const meshDirection<Type1,MeshType>&,                                       \
    const tmp<meshDirection<Type2,MeshType>>&                                   \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>&,                                  \
    const meshDirection<Type2,MeshType>&                                        \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>&,                                  \
    const tmp<meshDirection<Type2,MeshType>>&                                   \
);

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshDirection<ReturnType,MeshType>&,                                        \
    const Type1&,                                                               \
    const meshDirection<Type2,MeshType>&                                        \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const Type1&,                                                               \
    const meshDirection<Type2,MeshType>&                                        \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const Type1&,                                                               \
    const tmp<meshDirection<Type2,MeshType>>&                                   \
);

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshDirection<ReturnType,MeshType>&,                                        \
    const meshDirection<Type1,MeshType>&,                                       \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const meshDirection<Type1,MeshType>&,                                       \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshDirection<ReturnType,MeshType>> operator Op                             \
(                                                                               \
    const tmp<meshDirection<Type1,MeshType>>&,                                  \
    const Type2&                                                                \
);

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)               \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)
