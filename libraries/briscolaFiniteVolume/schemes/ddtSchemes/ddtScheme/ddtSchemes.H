#ifndef ddtSchemes_H
#define ddtSchemes_H

#include "ddtScheme.H"

#define makeDdtScheme(TYPE,MESHTYPE)                                            \
                                                                                \
    typedef ddtScheme<TYPE,MESHTYPE> ddtScheme##TYPE##MESHTYPE;                 \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        ddtScheme##TYPE##MESHTYPE,                                              \
        "ddtScheme",                                                            \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        ddtScheme##TYPE##MESHTYPE,                                              \
        dictionary                                                              \
    );


#define makeDdtSchemeType(SCHEME,TYPE,MESHTYPE)                                 \
                                                                                \
    typedef SCHEME##DdtScheme<TYPE,MESHTYPE>                                    \
        SCHEME##DdtScheme##TYPE##MESHTYPE;                                      \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        SCHEME##DdtScheme##TYPE##MESHTYPE,                                      \
        #SCHEME,                                                                \
        0                                                                       \
    );                                                                          \
                                                                                \
    ddtScheme<TYPE,MESHTYPE>::                                                  \
    adddictionaryConstructorToTable<SCHEME##DdtScheme##TYPE##MESHTYPE>          \
    add##SCHEME##DdtScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
