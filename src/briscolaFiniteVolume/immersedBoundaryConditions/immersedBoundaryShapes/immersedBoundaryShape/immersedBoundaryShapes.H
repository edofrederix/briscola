#ifndef immersedBoundaryShapes_H
#define immersedBoundaryShapes_H

#include "immersedBoundaryShape.H"

#define makeShapeType(SHAPE)                                                    \
                                                                                \
    defineTypeNameAndDebug                                                      \
    (                                                                           \
        SHAPE,                                                                  \
        0                                                                       \
    );                                                                          \
                                                                                \
    immersedBoundaryShape::adddictionaryConstructorToTable<SHAPE>               \
    add##SHAPE##ConstructorToTable_;

#endif