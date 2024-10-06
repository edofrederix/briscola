#ifndef prolongationSchemes_H
#define prolongationSchemes_H

#include "prolongationScheme.H"

#define makeProlongationScheme(TYPE,MESHTYPE)                                  \
                                                                               \
    typedef prolongationScheme<TYPE,MESHTYPE>                                  \
        prolongationScheme##TYPE##MESHTYPE;                                    \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        prolongationScheme##TYPE##MESHTYPE,                                    \
        "prolongationScheme",                                                  \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        prolongationScheme##TYPE##MESHTYPE,                                    \
        Istream                                                                \
    );

#define makeProlongationSchemeType(SCHEME,TYPE,MESHTYPE)                       \
                                                                               \
    typedef SCHEME##ProlongationScheme<TYPE,MESHTYPE>                          \
        SCHEME##ProlongationScheme##TYPE##MESHTYPE;                            \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SCHEME##ProlongationScheme##TYPE##MESHTYPE,                            \
        #SCHEME,                                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    prolongationScheme<TYPE,MESHTYPE>::                                        \
    addIstreamConstructorToTable                                               \
    <                                                                          \
        SCHEME##ProlongationScheme##TYPE##MESHTYPE                             \
    >                                                                          \
    add##SCHEME##ProlongationScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
