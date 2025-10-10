#ifndef faceDotGradientSchemes_H
#define faceDotGradientSchemes_H

#include "faceDotGradientScheme.H"

#define makeFaceDotGradientScheme(TYPE,MESHTYPE)                               \
                                                                               \
    typedef faceDotGradientScheme<TYPE,MESHTYPE>                               \
        faceDotGradientScheme##TYPE##MESHTYPE;                                 \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        faceDotGradientScheme##TYPE##MESHTYPE,                                 \
        "faceDotGradientScheme",                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        faceDotGradientScheme##TYPE##MESHTYPE,                                 \
        Istream                                                                \
    );


#define makeFaceDotGradientSchemeType(SCHEME,TYPE,MESHTYPE)                    \
                                                                               \
    typedef SCHEME##FaceDotGradientScheme<TYPE,MESHTYPE>                       \
        SCHEME##FaceDotGradientScheme##TYPE##MESHTYPE;                         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SCHEME##FaceDotGradientScheme##TYPE##MESHTYPE,                         \
        #SCHEME,                                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    faceDotGradientScheme<TYPE,MESHTYPE>::                                     \
    addIstreamConstructorToTable                                               \
    <                                                                          \
        SCHEME##FaceDotGradientScheme##TYPE##MESHTYPE                          \
    >                                                                          \
    add##SCHEME##FaceDotGradientScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
