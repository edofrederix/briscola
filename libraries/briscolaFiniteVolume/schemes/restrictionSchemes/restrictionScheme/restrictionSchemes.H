#ifndef restrictionSchemes_H
#define restrictionSchemes_H

#include "restrictionScheme.H"

#define makeRestrictionScheme(TYPE,MESHTYPE)                                    \
                                                                                \
    typedef restrictionScheme<TYPE,MESHTYPE>                                    \
        restrictionScheme##TYPE##MESHTYPE;                                      \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        restrictionScheme##TYPE##MESHTYPE,                                      \
        "restrictionScheme",                                                    \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        restrictionScheme##TYPE##MESHTYPE,                                      \
        dictionary                                                              \
    );

#define makeRestrictionSchemeType(SCHEME,TYPE,MESHTYPE)                         \
                                                                                \
    typedef SCHEME##RestrictionScheme<TYPE,MESHTYPE>                            \
        SCHEME##RestrictionScheme##TYPE##MESHTYPE;                              \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        SCHEME##RestrictionScheme##TYPE##MESHTYPE,                              \
        #SCHEME,                                                                \
        0                                                                       \
    );                                                                          \
                                                                                \
    restrictionScheme<TYPE,MESHTYPE>::                                          \
    adddictionaryConstructorToTable<SCHEME##RestrictionScheme##TYPE##MESHTYPE>  \
    add##SCHEME##RestrictionScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
