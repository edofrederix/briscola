#ifndef prolongationSchemes_H
#define prolongationSchemes_H

#include "prolongationScheme.H"

#define makeProlongationScheme(TYPE,MESHTYPE,DEFAULT)                          \
                                                                               \
    typedef prolongationScheme<TYPE,MESHTYPE>                                  \
        prolongationScheme##TYPE##MESHTYPE;                                    \
                                                                               \
    template<>                                                                 \
    const word prolongationScheme##TYPE##MESHTYPE::defaultScheme = DEFAULT;    \
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
