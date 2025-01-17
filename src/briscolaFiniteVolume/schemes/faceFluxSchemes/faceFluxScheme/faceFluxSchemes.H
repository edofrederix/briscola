#ifndef faceFluxSchemes_H
#define faceFluxSchemes_H

#include "faceFluxScheme.H"

#define makeFaceFluxSchemeType(SCHEME)                                         \
                                                                               \
    defineTypeNameAndDebug                                                     \
    (                                                                          \
        SCHEME##FaceFluxScheme,                                                \
        0                                                                      \
    );                                                                         \
                                                                               \
    faceFluxScheme::addIstreamConstructorToTable<SCHEME##FaceFluxScheme>       \
    add##SCHEME##FaceFluxScheme##ConstructorToTable_;

#endif
