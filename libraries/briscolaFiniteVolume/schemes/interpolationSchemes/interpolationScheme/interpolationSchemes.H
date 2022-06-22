#ifndef interpolationSchemes_H
#define interpolationSchemes_H

#include "interpolationScheme.H"

#define makeInterpolationScheme(TYPE,MESHTYPE)                                  \
                                                                                \
    typedef interpolationScheme<TYPE,MESHTYPE>                                  \
        interpolationScheme##TYPE##MESHTYPE;                                    \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        interpolationScheme##TYPE##MESHTYPE,                                    \
        "interpolationScheme",                                                  \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        interpolationScheme##TYPE##MESHTYPE,                                    \
        dictionary                                                              \
    );


#define makeInterpolationSchemeType(SCHEME,TYPE,MESHTYPE)                       \
                                                                                \
    typedef SCHEME##InterpolationScheme<TYPE,MESHTYPE>                          \
        SCHEME##InterpolationScheme##TYPE##MESHTYPE;                            \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        SCHEME##InterpolationScheme##TYPE##MESHTYPE,                            \
        #SCHEME,                                                                \
        0                                                                       \
    );                                                                          \
                                                                                \
    interpolationScheme<TYPE,MESHTYPE>::                                        \
    adddictionaryConstructorToTable                                             \
    <                                                                           \
        SCHEME##InterpolationScheme##TYPE##MESHTYPE                             \
    >                                                                           \
    add##SCHEME##InterpolationScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
