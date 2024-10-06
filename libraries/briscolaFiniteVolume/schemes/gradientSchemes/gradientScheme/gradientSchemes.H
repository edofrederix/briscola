#ifndef gradientSchemes_H
#define gradientSchemes_H

#include "gradientScheme.H"

#define makeGradientScheme(TYPE,MESHTYPE)                                       \
                                                                                \
    typedef gradientScheme<TYPE,MESHTYPE>                                       \
        gradientScheme##TYPE##MESHTYPE;                                         \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        gradientScheme##TYPE##MESHTYPE,                                         \
        "gradientScheme",                                                       \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        gradientScheme##TYPE##MESHTYPE,                                         \
        Istream                                                                 \
    );


#define makeGradientSchemeType(SCHEME,TYPE,MESHTYPE)                            \
                                                                                \
    typedef SCHEME##GradientScheme<TYPE,MESHTYPE>                               \
        SCHEME##GradientScheme##TYPE##MESHTYPE;                                 \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        SCHEME##GradientScheme##TYPE##MESHTYPE,                                 \
        #SCHEME,                                                                \
        0                                                                       \
    );                                                                          \
                                                                                \
    gradientScheme<TYPE,MESHTYPE>::                                             \
    addIstreamConstructorToTable<SCHEME##GradientScheme##TYPE##MESHTYPE>        \
    add##SCHEME##GradientScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
