#ifndef surfaceTensionSchemes_H
#define surfaceTensionSchemes_H

#include "surfaceTensionScheme.H"

#define makeSurfaceTensionScheme(SHORTNAME,SCHEME,SIGMA)                        \
                                                                                \
    typedef SCHEME<SIGMA> SCHEME##SIGMA;                                        \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        SCHEME##SIGMA,                                                          \
        #SHORTNAME,                                                             \
        0                                                                       \
    );                                                                          \
                                                                                \
    surfaceTensionScheme::adddictionaryConstructorToTable                       \
    <                                                                           \
        SCHEME##SIGMA                                                           \
    >                                                                           \
    SCHEME##SIGMA##ConstructorToTable_;

#endif