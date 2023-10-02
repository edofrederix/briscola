// Macros to iterative over levels, directions and internal cells

#define forAllLevels(FIELD, l, d, i, j, k)                                      \
    for (int l = 0; l < (FIELD).numberOfLevels(); l++)                          \
        forAllDirections((FIELD).level(l), d, i, j, k)

#define forAllLevelsReversed(FIELD, l, d, i, j, k)                              \
    for (int l = (FIELD).numberOfLevels()-1; l >= 0; l--)                       \
        forAllDirectionsReversed((FIELD).level(l), d, i, j, k)

#define UNARY_FUNCTION(ReturnType, Type1, Func)                                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshField<ReturnType,MeshType>&,                                            \
    const meshField<Type1,MeshType>&                                            \
);                                                                              \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>>                                             \
Func(const meshField<Type1,MeshType>&);                                         \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>>                                             \
Func(const tmp<meshField<Type1,MeshType>>&);

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshField<ReturnType,MeshType>&,                                            \
    const meshField<Type1,MeshType>&                                            \
);                                                                              \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshField<Type1,MeshType>&                                            \
);                                                                              \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshField<Type1,MeshType>>&                                       \
);

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                         \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshField<ReturnType,MeshType>&,                                            \
    const meshField<Type1,MeshType>&,                                           \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshField<Type1,MeshType>&,                                           \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshField<Type1,MeshType>&,                                           \
    const tmp<meshField<Type2,MeshType>>&                                       \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshField<Type1,MeshType>>&,                                      \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshField<Type1,MeshType>>&,                                      \
    const tmp<meshField<Type2,MeshType>>&                                       \
);

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshField<ReturnType,MeshType>&,                                            \
    const Type1&,                                                               \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const Type1&,                                                               \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const Type1&,                                                               \
    const tmp<meshField<Type2,MeshType>>&                                       \
);

#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                 \
                                                                                \
TEMPLATE                                                                        \
void Func                                                                       \
(                                                                               \
    meshField<ReturnType,MeshType>&,                                            \
    const meshField<Type1,MeshType>&,                                           \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const meshField<Type1,MeshType>&,                                           \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> Func                                        \
(                                                                               \
    const tmp<meshField<Type1,MeshType>>&,                                      \
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
    meshField<ReturnType,MeshType>&,                                            \
    const meshField<Type1,MeshType>&,                                           \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshField<Type1,MeshType>&,                                           \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshField<Type1,MeshType>&,                                           \
    const tmp<meshField<Type2,MeshType>>&                                       \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshField<Type1,MeshType>>&,                                      \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshField<Type1,MeshType>>&,                                      \
    const tmp<meshField<Type2,MeshType>>&                                       \
);

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshField<ReturnType,MeshType>&,                                            \
    const Type1&,                                                               \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const Type1&,                                                               \
    const meshField<Type2,MeshType>&                                            \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const Type1&,                                                               \
    const tmp<meshField<Type2,MeshType>>&                                       \
);

#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)           \
                                                                                \
TEMPLATE                                                                        \
void OpFunc                                                                     \
(                                                                               \
    meshField<ReturnType,MeshType>&,                                            \
    const meshField<Type1,MeshType>&,                                           \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const meshField<Type1,MeshType>&,                                           \
    const Type2&                                                                \
);                                                                              \
                                                                                \
TEMPLATE                                                                        \
tmp<meshField<ReturnType,MeshType>> operator Op                                 \
(                                                                               \
    const tmp<meshField<Type1,MeshType>>&,                                      \
    const Type2&                                                                \
);

#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)               \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)
