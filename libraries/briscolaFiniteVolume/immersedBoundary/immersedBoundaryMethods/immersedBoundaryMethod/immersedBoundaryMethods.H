#ifndef immersedBoundaryMethods_H
#define immersedBoundaryMethods_H

#include "immersedBoundaryMethod.H"

#define makeIBM(TYPE,MESHTYPE)                                                  \
                                                                                \
    typedef immersedBoundaryMethod<TYPE,MESHTYPE>                               \
        immersedBoundaryMethod##TYPE##MESHTYPE;                                 \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        immersedBoundaryMethod##TYPE##MESHTYPE,                                 \
        "immersedBoundaryMethod",                                               \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        immersedBoundaryMethod##TYPE##MESHTYPE,                                 \
        dictionary                                                              \
    );


#define makeIBMType(IBM,TYPE,MESHTYPE)                                          \
                                                                                \
    typedef IBM<TYPE,MESHTYPE>                                                  \
        IBM##TYPE##MESHTYPE;                                                    \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        IBM##TYPE##MESHTYPE,                                                    \
        #IBM,                                                                   \
        0                                                                       \
    );                                                                          \
                                                                                \
    immersedBoundaryMethod<TYPE,MESHTYPE>::                                     \
        adddictionaryConstructorToTable<IBM##TYPE##MESHTYPE>                    \
            add##IBM##TYPE##MESHTYPE##ConstructorToTable_;

#endif
