#define UNARY_FUNCTION(ReturnType, Type1, Func)                                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const meshLevel<Type1,MeshType>&                                            \
);                                                                              \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>>                                             \
Func(const meshLevel<Type1,MeshType>&);                                         \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>>                                             \
Func(const tmp<meshLevel<Type1,MeshType>>&);

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const meshLevel<Type1,MeshType>&                                            \
);                                                                              \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>>                                             \
operator Op(const meshLevel<Type1,MeshType>&);                                  \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>>                                             \
operator Op(const tmp<meshLevel<Type1,MeshType>>&);

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                         \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const meshLevel<Type1,MeshType>&,                                           \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshLevel<Type1,MeshType>&,                                           \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshLevel<Type1,MeshType>&,                                           \
    const tmp<meshLevel<Type2,MeshType>>&                                       \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>&,                                      \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>&,                                      \
    const tmp<meshLevel<Type2,MeshType>>&                                       \
);

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const Type1&,                                                               \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const Type1&,                                                               \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const Type1&,                                                               \
    const tmp<meshLevel<Type2,MeshType>>&                                       \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const List<Type1>&,                                                         \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const List<Type1>&,                                                         \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const List<Type1>&,                                                         \
    const tmp<meshLevel<Type2,MeshType>>&                                       \
);

#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const meshLevel<Type1,MeshType>&,                                           \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshLevel<Type1,MeshType>&,                                           \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>&,                                      \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const meshLevel<Type1,MeshType>&,                                           \
    const List<Type2>&                                                          \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshLevel<Type1,MeshType>&,                                           \
    const List<Type2>&                                                          \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>&,                                      \
    const List<Type2>&                                                          \
);

#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                     \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                   \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const meshLevel<Type1,MeshType>&,                                           \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshLevel<Type1,MeshType>&,                                           \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshLevel<Type1,MeshType>&,                                           \
    const tmp<meshLevel<Type2,MeshType>>&                                       \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>&,                                      \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>&,                                      \
    const tmp<meshLevel<Type2,MeshType>>&                                       \
);

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const Type1&,                                                               \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const Type1&,                                                               \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const Type1&,                                                               \
    const tmp<meshLevel<Type2,MeshType>>&                                       \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const List<Type1>&,                                                         \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const List<Type1>&,                                                         \
    const meshLevel<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const List<Type1>&,                                                         \
    const tmp<meshLevel<Type2,MeshType>>&                                       \
);

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const meshLevel<Type1,MeshType>&,                                           \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshLevel<Type1,MeshType>&,                                           \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>&,                                      \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshLevel<ReturnType,MeshType>&,                                            \
    const meshLevel<Type1,MeshType>&,                                           \
    const List<Type2>&                                                          \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshLevel<Type1,MeshType>&,                                           \
    const List<Type2>&                                                          \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshLevel<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshLevel<Type1,MeshType>>&,                                      \
    const List<Type2>&                                                          \
);

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)               \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)
