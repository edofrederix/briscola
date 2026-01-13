#ifndef restrictionSchemes_H
#define restrictionSchemes_H

#include "restrictionScheme.H"

// Make the base restriction scheme

#define makeRestrictionSchemeBase(TYPE,MESHTYPE,DEFAULT)                       \
                                                                               \
    typedef restrictionScheme<TYPE,MESHTYPE>                                   \
        restrictionScheme##TYPE##MESHTYPE;                                     \
                                                                               \
    template<>                                                                 \
    const word restrictionScheme##TYPE##MESHTYPE::defaultScheme = DEFAULT;     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        restrictionScheme##TYPE##MESHTYPE,                                     \
        "restrictionScheme",                                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        restrictionScheme##TYPE##MESHTYPE,                                     \
        Istream                                                                \
    );


// Make a restriction scheme templated for a type and mesh type

#define makeRestrictionScheme(SCHEME,TYPE,MESHTYPE)                            \
                                                                               \
    typedef SCHEME##RestrictionScheme<TYPE,MESHTYPE>                           \
        SCHEME##RestrictionScheme##TYPE##MESHTYPE;                             \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SCHEME##RestrictionScheme##TYPE##MESHTYPE,                             \
        #SCHEME,                                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    restrictionScheme<TYPE,MESHTYPE>::                                         \
    addIstreamConstructorToTable<SCHEME##RestrictionScheme##TYPE##MESHTYPE>    \
    add##SCHEME##RestrictionScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
