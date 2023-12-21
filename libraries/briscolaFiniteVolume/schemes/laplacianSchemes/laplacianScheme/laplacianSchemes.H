#ifndef laplacianSchemes_H
#define laplacianSchemes_H

#include "laplacianScheme.H"

#define makeLaplacianScheme(STYPE,TYPE,MESHTYPE)                               \
                                                                               \
    typedef laplacianScheme<STYPE,TYPE,MESHTYPE>                               \
        laplacianScheme##STYPE##TYPE##MESHTYPE;                                \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        laplacianScheme##STYPE##TYPE##MESHTYPE,                                \
        "laplacianScheme",                                                     \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        laplacianScheme##STYPE##TYPE##MESHTYPE,                                \
        dictionary                                                             \
    );


#define makeLaplacianSchemeType(SCHEME,STYPE,TYPE,MESHTYPE)                    \
                                                                               \
    typedef SCHEME##LaplacianScheme<STYPE,TYPE,MESHTYPE>                       \
        SCHEME##LaplacianScheme##STYPE##TYPE##MESHTYPE;                        \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SCHEME##LaplacianScheme##STYPE##TYPE##MESHTYPE,                        \
        #SCHEME,                                                               \
        0                                                                      \
    );

#define addSpecificLaplacianSchemeType(SCHEME,STYPE,TYPE,MESHTYPE)             \
                                                                               \
    typedef SCHEME##LaplacianScheme<TYPE,MESHTYPE>                             \
        SCHEME##LaplacianScheme##STYPE##TYPE##MESHTYPE;                        \
                                                                               \
    laplacianScheme<STYPE,TYPE,MESHTYPE>::                                     \
    adddictionaryConstructorToTable                                            \
    <                                                                          \
        SCHEME##LaplacianScheme##STYPE##TYPE##MESHTYPE                         \
    >                                                                          \
    add##SCHEME##LaplacianScheme##STYPE##TYPE##MESHTYPE##ConstructorToTable_;

#endif
