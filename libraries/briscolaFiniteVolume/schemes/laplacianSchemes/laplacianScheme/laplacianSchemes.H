#ifndef laplacianSchemes_H
#define laplacianSchemes_H

#include "laplacianScheme.H"

#define makeLaplacianScheme(TYPE,MESHTYPE)                                      \
                                                                                \
    typedef laplacianScheme<TYPE,MESHTYPE>                                      \
        laplacianScheme##TYPE##MESHTYPE;                                        \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        laplacianScheme##TYPE##MESHTYPE,                                        \
        "laplacianScheme",                                                      \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        laplacianScheme##TYPE##MESHTYPE,                                        \
        dictionary                                                              \
    );


#define makeLaplacianSchemeType(SCHEME,TYPE,MESHTYPE)                           \
                                                                                \
    typedef SCHEME##LaplacianScheme<TYPE,MESHTYPE>                              \
        SCHEME##LaplacianScheme##TYPE##MESHTYPE;                                \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        SCHEME##LaplacianScheme##TYPE##MESHTYPE,                                \
        #SCHEME,                                                                \
        0                                                                       \
    );                                                                          \
                                                                                \
    laplacianScheme<TYPE,MESHTYPE>::                                            \
    adddictionaryConstructorToTable<SCHEME##LaplacianScheme##TYPE##MESHTYPE>    \
    add##SCHEME##LaplacianScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
