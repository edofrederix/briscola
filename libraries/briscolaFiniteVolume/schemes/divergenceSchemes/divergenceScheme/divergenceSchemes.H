#ifndef divergenceSchemes_H
#define divergenceSchemes_H

#include "divergenceScheme.H"

#define makeDivergenceScheme(STYPE,TYPE,MESHTYPE)                              \
                                                                               \
    typedef divergenceScheme<STYPE,TYPE,MESHTYPE>                              \
        divergenceScheme##STYPE##TYPE##MESHTYPE;                               \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        divergenceScheme##STYPE##TYPE##MESHTYPE,                               \
        "divergenceScheme",                                                    \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        divergenceScheme##STYPE##TYPE##MESHTYPE,                               \
        Istream                                                                \
    );

#define makeDivergenceSchemeType(SCHEME,STYPE,TYPE,MESHTYPE)                   \
                                                                               \
    typedef SCHEME##DivergenceScheme<TYPE,MESHTYPE>                            \
        SCHEME##DivergenceScheme##STYPE##TYPE##MESHTYPE;                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SCHEME##DivergenceScheme##STYPE##TYPE##MESHTYPE,                       \
        #SCHEME,                                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    divergenceScheme<STYPE,TYPE,MESHTYPE>::                                    \
    addIstreamConstructorToTable                                               \
    <                                                                          \
        SCHEME##DivergenceScheme##STYPE##TYPE##MESHTYPE                        \
    >                                                                          \
    add##SCHEME##DivergenceScheme##STYPE##TYPE##MESHTYPE##ConstructorToTable_;

#endif
