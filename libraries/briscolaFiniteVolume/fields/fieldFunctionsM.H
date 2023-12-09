#ifndef fieldFunctionsM_H
#define fieldFunctionsM_H

// Specialized iterators are defined for directions, levels and fields. These
// are then combined in an overloaded iterator (forAllCells for cells and
// forAllFaces for faces), selecting the specialized one based on the number of
// arguments.

// For all cells

#define forAllCellsInDirection(DIR,i,j,k)                                      \
    for (int i = (DIR).I().left();   i < (DIR).I().right(); i++)               \
    for (int j = (DIR).I().bottom(); j < (DIR).I().top();   j++)               \
    for (int k = (DIR).I().aft();    k < (DIR).I().fore();  k++)

#define forAllCellsInLevel(LEV,d,i,j,k)                                        \
    for (int d = 0; d < (LEV).numberOfDirections(); d++)                       \
        forAllCellsInDirection((LEV).direction(d), i, j, k)

#define forAllCellsInField(FIELD,l,d,i,j,k)                                    \
    for (int l = 0; l < (FIELD).numberOfLevels(); l++)                         \
        forAllCellsInLevel((FIELD).level(l), d, i, j, k)

#define GET_CELL_ITERATOR(_1,_2,_3,_4,_5,_6,NAME,...) NAME
#define forAllCells(...)                                                       \
    GET_CELL_ITERATOR                                                          \
    (                                                                          \
        __VA_ARGS__,                                                           \
        forAllCellsInField,                                                    \
        forAllCellsInLevel,                                                    \
        forAllCellsInDirection                                                 \
    )(__VA_ARGS__)

// For all faces

#define forAllFacesInDirection(DIR,fd,i,j,k)                                   \
    for (int fd = 0; fd < 3; fd++)                                             \
    for (int i = (DIR).I().left();   i < (DIR).I().right() + (fd == 0); i++)   \
    for (int j = (DIR).I().bottom(); j < (DIR).I().top()   + (fd == 1); j++)   \
    for (int k = (DIR).I().aft();    k < (DIR).I().fore()  + (fd == 2); k++)

#define forAllFacesInLevel(LEV,d,fd,i,j,k)                                     \
    for (int d = 0; d < (LEV).numberOfDirections(); d++)                       \
        forAllFacesInDirection((LEV).direction(d), fd, i, j, k)

#define forAllFacesInField(FIELD,l,d,fd,i,j,k)                                 \
    for (int l = 0; l < (FIELD).numberOfLevels(); l++)                         \
        forAllFacesInLevel((FIELD).level(l), d, fd, i, j, k)

#define GET_FACE_ITERATOR(_1,_2,_3,_4,_5,_6,_7,NAME,...) NAME
#define forAllFaces(...)                                                       \
    GET_FACE_ITERATOR                                                          \
    (                                                                          \
        __VA_ARGS__,                                                           \
        forAllFacesInField,                                                    \
        forAllFacesInLevel,                                                    \
        forAllFacesInDirection                                                 \
    )(__VA_ARGS__)

#endif
