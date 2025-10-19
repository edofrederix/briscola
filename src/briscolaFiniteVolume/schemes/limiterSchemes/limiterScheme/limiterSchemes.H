#ifndef limiterSchemes_H
#define limiterSchemes_H

#include "limiterScheme.H"

#define makeLimiterScheme(TYPE,MESHTYPE)                                       \
                                                                               \
    typedef limiterScheme<TYPE,MESHTYPE>                                       \
        limiterScheme##TYPE##MESHTYPE;                                         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        limiterScheme##TYPE##MESHTYPE,                                         \
        "limiterScheme",                                                       \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        limiterScheme##TYPE##MESHTYPE,                                         \
        Istream                                                                \
    );                                                                         \
                                                                               \
/* Instantiate the r-function */                                               \
                                                                               \
template<>                                                                     \
tmp<faceField<scalar,MESHTYPE>>                                                \
limiterScheme<TYPE,MESHTYPE>::r                                                \
(                                                                              \
    const faceField<scalar,MESHTYPE>& phi,                                     \
    const meshField<TYPE,MESHTYPE>& field                                      \
)                                                                              \
{                                                                              \
    return limiterScheme<TYPE,MESHTYPE>::rType(phi,field);                     \
}

#define makeLimiterSchemeType(SCHEME,TYPE,MESHTYPE)                            \
                                                                               \
    typedef SCHEME##LimiterScheme<TYPE,MESHTYPE>                               \
        SCHEME##LimiterScheme##TYPE##MESHTYPE;                                 \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SCHEME##LimiterScheme##TYPE##MESHTYPE,                                 \
        #SCHEME,                                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    limiterScheme<TYPE,MESHTYPE>::                                             \
    addIstreamConstructorToTable                                               \
    <                                                                          \
        SCHEME##LimiterScheme##TYPE##MESHTYPE                                  \
    >                                                                          \
    add##SCHEME##LimiterScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
