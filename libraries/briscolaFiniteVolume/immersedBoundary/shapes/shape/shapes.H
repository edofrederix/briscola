#ifndef shapes_H
#define shapes_H

#include "shape.H"

#define makeShapeType(SHAPE)                                                    \
                                                                                \
    defineTypeNameAndDebug                                                      \
    (                                                                           \
        SHAPE,                                                                  \
        0                                                                       \
    );                                                                          \
                                                                                \
    shape::adddictionaryConstructorToTable<SHAPE>                               \
    add##SHAPE##ConstructorToTable_;

#endif