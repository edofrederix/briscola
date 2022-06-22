#ifndef reconstructionSchemes_H
#define reconstructionSchemes_H

#include "reconstructionScheme.H"

#define makeReconstructionScheme(TYPE)                                          \
                                                                                \
    typedef reconstructionScheme<TYPE>                                          \
        reconstructionScheme##TYPE;                                             \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        reconstructionScheme##TYPE,                                             \
        "reconstructionScheme",                                                 \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        reconstructionScheme##TYPE,                                             \
        dictionary                                                              \
    );


#define makeReconstructionSchemeType(SCHEME,TYPE)                               \
                                                                                \
    typedef SCHEME##ReconstructionScheme<TYPE>                                  \
        SCHEME##ReconstructionScheme##TYPE;                                     \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        SCHEME##ReconstructionScheme##TYPE,                                     \
        #SCHEME,                                                                \
        0                                                                       \
    );                                                                          \
                                                                                \
    reconstructionScheme<TYPE>::                                                \
    adddictionaryConstructorToTable<SCHEME##ReconstructionScheme##TYPE>         \
    add##SCHEME##ReconstructionScheme##TYPE##ConstructorToTable_;

#endif
