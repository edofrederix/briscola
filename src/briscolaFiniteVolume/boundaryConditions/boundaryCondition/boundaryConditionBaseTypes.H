#ifndef boundaryConditionBaseTypes_H
#define boundaryConditionBaseTypes_H

#include "NamedEnum.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Enumeration listing the base boundary condition types
enum boundaryConditionBaseType
{
    DUMMYBC,
    EMPTYBC,
    PARALLELBC,
    PERIODICBC,
    DIRICHLETBC,
    NEUMANNBC,
    ROBINBC
};

const NamedEnum<boundaryConditionBaseType,7> boundaryConditionBaseTypeNames;

}

}

}

#endif
