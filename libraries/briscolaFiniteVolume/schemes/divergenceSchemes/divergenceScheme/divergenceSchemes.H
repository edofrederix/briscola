#ifndef divergenceSchemes_H
#define divergenceSchemes_H

#include "divergenceScheme.H"

#define makeDivergenceScheme(TYPE,MESHTYPE)                                    \
                                                                               \
    typedef divergenceScheme<TYPE,MESHTYPE>                                    \
        divergenceScheme##TYPE##MESHTYPE;                                      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        divergenceScheme##TYPE##MESHTYPE,                                      \
        "divergenceScheme",                                                    \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        divergenceScheme##TYPE##MESHTYPE,                                      \
        dictionary                                                             \
    );

#define makeDivergenceSchemeType(SCHEME,TYPE,MESHTYPE)                         \
                                                                               \
    typedef SCHEME##DivergenceScheme<TYPE,MESHTYPE>                            \
        SCHEME##DivergenceScheme##TYPE##MESHTYPE;                              \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SCHEME##DivergenceScheme##TYPE##MESHTYPE,                              \
        #SCHEME,                                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    divergenceScheme<TYPE,MESHTYPE>::                                          \
    adddictionaryConstructorToTable<SCHEME##DivergenceScheme##TYPE##MESHTYPE>  \
    add##SCHEME##DivergenceScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
