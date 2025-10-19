#ifndef faceRestrictionSchemes_H
#define faceRestrictionSchemes_H

#include "faceRestrictionScheme.H"

// Make the base restriction scheme

#define makeFaceRestrictionSchemeBase(TYPE,MESHTYPE,DEFAULT)                   \
                                                                               \
    typedef faceRestrictionScheme<TYPE,MESHTYPE>                               \
        faceRestrictionScheme##TYPE##MESHTYPE;                                 \
                                                                               \
    template<>                                                                 \
    const word faceRestrictionScheme##TYPE##MESHTYPE::defaultScheme = DEFAULT; \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        faceRestrictionScheme##TYPE##MESHTYPE,                                 \
        "faceRestrictionScheme",                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        faceRestrictionScheme##TYPE##MESHTYPE,                                 \
        Istream                                                                \
    );


// Make a restriction scheme templated for a type and mesh type

#define makeFaceRestrictionScheme(SCHEME,TYPE,MESHTYPE)                        \
                                                                               \
    typedef SCHEME##FaceRestrictionScheme<TYPE,MESHTYPE>                       \
        SCHEME##FaceRestrictionScheme##TYPE##MESHTYPE;                         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SCHEME##FaceRestrictionScheme##TYPE##MESHTYPE,                         \
        #SCHEME,                                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    faceRestrictionScheme<TYPE,MESHTYPE>::                                     \
    addIstreamConstructorToTable<SCHEME##FaceRestrictionScheme##TYPE##MESHTYPE>\
    add##SCHEME##FaceRestrictionScheme##TYPE##MESHTYPE##ConstructorToTable_;

#endif
