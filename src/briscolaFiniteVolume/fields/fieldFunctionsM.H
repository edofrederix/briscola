#ifndef fieldFunctionsM_H
#define fieldFunctionsM_H

// Specialized iterators are defined for directions, levels and fields. These
// are then combined in an overloaded iterator (forAllCells for cells and
// forAllFaces for faces), selecting the specialized one based on the number of
// arguments.

// For all cells

#define forAllCellsInDirection(DIR,i,j,k)                                      \
    for (int i = (DIR).I()[0]; i < (DIR).I()[1]; i++)                          \
    for (int j = (DIR).I()[2]; j < (DIR).I()[3]; j++)                          \
    for (int k = (DIR).I()[4]; k < (DIR).I()[5]; k++)

#define forAllCellsInLevel(LEV,d,i,j,k)                                        \
    for (int d = 0; d < (LEV).numberOfDirections(); d++)                       \
        forAllCellsInDirection((LEV).direction(d), i, j, k)

#define forAllCellsInField(FIE,l,d,i,j,k)                                      \
    for (int l = 0; l < (FIE).numberOfLevels(); l++)                           \
        forAllCellsInLevel((FIE).level(l), d, i, j, k)

#define GET_CELL_ITERATOR(_1,_2,_3,_4,_5,_6,NAME,...) NAME
#define forAllCells(...)                                                       \
    GET_CELL_ITERATOR                                                          \
    (                                                                          \
        __VA_ARGS__,                                                           \
        forAllCellsInField,                                                    \
        forAllCellsInLevel,                                                    \
        forAllCellsInDirection                                                 \
    )(__VA_ARGS__)

#define forAllCellsReversedInDirection(DIR,i,j,k)                              \
    for (int i = (DIR).I()[1] - 1; i >= (DIR).I()[0]; i--)                     \
    for (int j = (DIR).I()[3] - 1; j >= (DIR).I()[2]; j--)                     \
    for (int k = (DIR).I()[5] - 1; k >= (DIR).I()[4]; k--)

#define forAllCellsReversedInLevel(LEV,d,i,j,k)                                \
    for (int d = 0; d < (LEV).numberOfDirections(); d++)                       \
        forAllCellsReversedInDirection((LEV).direction(d), i, j, k)

#define forAllCellsReversedInField(FIE,l,d,i,j,k)                              \
    for (int l = 0; l < (FIE).numberOfLevels(); l++)                           \
        forAllCellsReversedInLevel((FIE).level(l), d, i, j, k)

#define GET_CELL_REV_ITERATOR(_1,_2,_3,_4,_5,_6,NAME,...) NAME
#define forAllCellsReversed(...)                                               \
    GET_CELL_REV_ITERATOR                                                      \
    (                                                                          \
        __VA_ARGS__,                                                           \
        forAllCellsReversedInField,                                            \
        forAllCellsReversedInLevel,                                            \
        forAllCellsReversedInDirection                                         \
    )(__VA_ARGS__)

// For all faces (takes only a meshFieldSet as argument)

#define forAllFacesInDirection(SET,fd,i,j,k)                                   \
    forAll((SET), fd)                                                          \
    for (int i=(SET)[fd][0][0].I()[0]; i<(SET)[fd][0][0].I()[1]+(fd==0); i++)  \
    for (int j=(SET)[fd][0][0].I()[2]; j<(SET)[fd][0][0].I()[3]+(fd==1); j++)  \
    for (int k=(SET)[fd][0][0].I()[4]; k<(SET)[fd][0][0].I()[5]+(fd==2); k++)

#define forAllFacesInLevel(SET,fd,d,i,j,k)                                     \
    forAll((SET), fd)                                                          \
    forAll((SET)[fd][0], d)                                                    \
    for (int i=(SET)[fd][0][d].I()[0]; i<(SET)[fd][0][d].I()[1]+(fd==0); i++)  \
    for (int j=(SET)[fd][0][d].I()[2]; j<(SET)[fd][0][d].I()[3]+(fd==1); j++)  \
    for (int k=(SET)[fd][0][d].I()[4]; k<(SET)[fd][0][d].I()[5]+(fd==2); k++)

#define forAllFacesInField(SET,fd,l,d,i,j,k)                                   \
    forAll((SET), fd)                                                          \
    forAll((SET)[fd], l)                                                       \
    forAll((SET)[fd][l], d)                                                    \
    for (int i=(SET)[fd][l][d].I()[0]; i<(SET)[fd][l][d].I()[1]+(fd==0); i++)  \
    for (int j=(SET)[fd][l][d].I()[2]; j<(SET)[fd][l][d].I()[3]+(fd==1); j++)  \
    for (int k=(SET)[fd][l][d].I()[4]; k<(SET)[fd][l][d].I()[5]+(fd==2); k++)

#define GET_FACE_ITERATOR(_1,_2,_3,_4,_5,_6,_7,NAME,...) NAME
#define forAllFaces(...)                                                       \
    GET_FACE_ITERATOR                                                          \
    (                                                                          \
        __VA_ARGS__,                                                           \
        forAllFacesInField,                                                    \
        forAllFacesInLevel,                                                    \
        forAllFacesInDirection                                                 \
    )(__VA_ARGS__)

// For all specific faces (the face direction is set)

#define forAllSpecificFacesInDirection(SET,fd,i,j,k)                           \
    for (int i=(SET)[fd][0][0].I()[0]; i<(SET)[fd][0][0].I()[1]+(fd==0); i++)  \
    for (int j=(SET)[fd][0][0].I()[2]; j<(SET)[fd][0][0].I()[3]+(fd==1); j++)  \
    for (int k=(SET)[fd][0][0].I()[4]; k<(SET)[fd][0][0].I()[5]+(fd==2); k++)

#define forAllSpecificFacesInLevel(SET,fd,d,i,j,k)                             \
    forAll((SET)[fd][0], d)                                                    \
    for (int i=(SET)[fd][0][d].I()[0]; i<(SET)[fd][0][d].I()[1]+(fd==0); i++)  \
    for (int j=(SET)[fd][0][d].I()[2]; j<(SET)[fd][0][d].I()[3]+(fd==1); j++)  \
    for (int k=(SET)[fd][0][d].I()[4]; k<(SET)[fd][0][d].I()[5]+(fd==2); k++)

#define forAllSpecificFacesInField(SET,fd,l,d,i,j,k)                           \
    forAll((SET)[fd], l)                                                       \
    forAll((SET)[fd][l], d)                                                    \
    for (int i=(SET)[fd][l][d].I()[0]; i<(SET)[fd][l][d].I()[1]+(fd==0); i++)  \
    for (int j=(SET)[fd][l][d].I()[2]; j<(SET)[fd][l][d].I()[3]+(fd==1); j++)  \
    for (int k=(SET)[fd][l][d].I()[4]; k<(SET)[fd][l][d].I()[5]+(fd==2); k++)

#define GET_FACE_ITERATOR(_1,_2,_3,_4,_5,_6,_7,NAME,...) NAME
#define forAllSpecificFaces(...)                                               \
    GET_FACE_ITERATOR                                                          \
    (                                                                          \
        __VA_ARGS__,                                                           \
        forAllSpecificFacesInField,                                            \
        forAllSpecificFacesInLevel,                                            \
        forAllSpecificFacesInDirection                                         \
    )(__VA_ARGS__)

// For all faces in a specific level

#define forAllFacesInSpecificLevel(SET,fd,l,d,i,j,k)                           \
    forAll((SET), fd)                                                          \
    forAll((SET)[fd][l], d)                                                    \
    for (int i=(SET)[fd][l][d].I()[0]; i<(SET)[fd][l][d].I()[1]+(fd==0); i++)  \
    for (int j=(SET)[fd][l][d].I()[2]; j<(SET)[fd][l][d].I()[3]+(fd==1); j++)  \
    for (int k=(SET)[fd][l][d].I()[4]; k<(SET)[fd][l][d].I()[5]+(fd==2); k++)

// For all faces in a specific level and direction

#define forAllFacesInSpecificDirection(SET,fd,l,d,i,j,k)                       \
    forAll((SET), fd)                                                          \
    for (int i=(SET)[fd][l][d].I()[0]; i<(SET)[fd][l][d].I()[1]+(fd==0); i++)  \
    for (int j=(SET)[fd][l][d].I()[2]; j<(SET)[fd][l][d].I()[3]+(fd==1); j++)  \
    for (int k=(SET)[fd][l][d].I()[4]; k<(SET)[fd][l][d].I()[5]+(fd==2); k++)

// For all specific faces in a specific level

#define forAllSpecificFacesInSpecificLevel(SET,fd,l,d,i,j,k)                   \
    forAll((SET)[fd][l], d)                                                    \
    for (int i=(SET)[fd][l][d].I()[0]; i<(SET)[fd][l][d].I()[1]+(fd==0); i++)  \
    for (int j=(SET)[fd][l][d].I()[2]; j<(SET)[fd][l][d].I()[3]+(fd==1); j++)  \
    for (int k=(SET)[fd][l][d].I()[4]; k<(SET)[fd][l][d].I()[5]+(fd==2); k++)

// For all specific faces in a specific level and direction

#define forAllSpecificFacesInSpecificDirection(SET,fd,l,d,i,j,k)               \
    for (int i=(SET)[fd][l][d].I()[0]; i<(SET)[fd][l][d].I()[1]+(fd==0); i++)  \
    for (int j=(SET)[fd][l][d].I()[2]; j<(SET)[fd][l][d].I()[3]+(fd==1); j++)  \
    for (int k=(SET)[fd][l][d].I()[4]; k<(SET)[fd][l][d].I()[5]+(fd==2); k++)

#endif
