#ifndef faceGradientSchemes_H
#define faceGradientSchemes_H

#include "faceGradientScheme.H"

#define makeFaceGradientScheme(TYPE,MESHTYPE)                                   \
                                                                                \
    typedef faceGradientScheme<TYPE,MESHTYPE>                                   \
        faceGradientScheme##TYPE##MESHTYPE;                                     \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        faceGradientScheme##TYPE##MESHTYPE,                                     \
        "faceGradientScheme",                                                   \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        faceGradientScheme##TYPE##MESHTYPE,                                     \
        Istream                                                                 \
    );


#define makeFaceGradientSchemeType(SCHEME,TYPE,MESHTYPE)                        \
                                                                                \
    typedef SCHEME##FaceGradientScheme<TYPE,MESHTYPE>                           \
        SCHEME##FaceGradientScheme##TYPE##MESHTYPE;                             \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        SCHEME##FaceGradientScheme##TYPE##MESHTYPE,                             \
        #SCHEME,                                                                \
        0                                                                       \
    );                                                                          \
                                                                                \
    faceGradientScheme<TYPE,MESHTYPE>::                                         \
    addIstreamConstructorToTable<SCHEME##FaceGradientScheme##TYPE##MESHTYPE>    \
    add##SCHEME##FaceGradientScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
