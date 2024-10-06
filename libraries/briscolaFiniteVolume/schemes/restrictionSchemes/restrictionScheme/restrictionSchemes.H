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

// Make a restriction scheme templated for a mesh type

#define makeRestrictionSchemeMeshType(SCHEME,NAME,TYPE,MESHTYPE)               \
                                                                               \
    typedef SCHEME##RestrictionScheme<MESHTYPE>                                \
        SCHEME##RestrictionScheme##MESHTYPE;                                   \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SCHEME##RestrictionScheme##MESHTYPE,                                   \
        #NAME,                                                                 \
        0                                                                      \
    );                                                                         \
                                                                               \
    restrictionScheme<TYPE,MESHTYPE>::                                         \
    addIstreamConstructorToTable<SCHEME##RestrictionScheme##MESHTYPE>          \
    add##SCHEME##RestrictionScheme##MESHTYPE##ConstructorToTable_;

// Make a restriction scheme templated for a type

#define makeRestrictionSchemeType(SCHEME,NAME,TYPE,MESHTYPE)                   \
                                                                               \
    typedef SCHEME##RestrictionScheme<TYPE>                                    \
        SCHEME##RestrictionScheme##TYPE;                                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SCHEME##RestrictionScheme##TYPE,                                       \
        #NAME,                                                                 \
        0                                                                      \
    );                                                                         \
                                                                               \
    restrictionScheme<TYPE,MESHTYPE>::                                         \
    addIstreamConstructorToTable<SCHEME##RestrictionScheme##TYPE>              \
    add##SCHEME##RestrictionScheme##TYPE##ConstructorToTable_;

// Make an non-templated restriction scheme

#define makeRestrictionSchemeNoTemplate(SCHEME,TYPE,MESHTYPE)                  \
                                                                               \
    defineTypeNameAndDebug                                                     \
    (                                                                          \
        SCHEME##RestrictionScheme,                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    restrictionScheme<TYPE,MESHTYPE>::                                         \
    addIstreamConstructorToTable<SCHEME##RestrictionScheme>                    \
    add##SCHEME##RestrictionScheme##ConstructorToTable_;

#endif
