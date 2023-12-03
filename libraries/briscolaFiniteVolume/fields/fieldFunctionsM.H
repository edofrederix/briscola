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

#endif
