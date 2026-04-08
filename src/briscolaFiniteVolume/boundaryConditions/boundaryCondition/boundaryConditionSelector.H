#ifndef boundaryConditionSelector_H
#define boundaryConditionSelector_H

#include "boundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Select all boundary conditions

class allBcs
{
public:
    template<class Type, class MeshType>
    static bool match(const boundaryCondition<Type,MeshType>&)
    {
        return true;
    }
};

// Select boundary conditions with a certain boundary type

template<class Boundary>
class bcsOfType
{
public:
    template<class Type, class MeshType>
    static bool match(const boundaryCondition<Type,MeshType>& bc)
    {
        return bc.b().template castable<Boundary>();
    }
};

// Select eliminated boundary conditions

class eliminatedBcs
{
public:
    template<class Type, class MeshType>
    static bool match(const boundaryCondition<Type,MeshType>& bc)
    {
        return bc.eliminated();
    }
};

// Select non-eliminated boundary conditions

class nonEliminatedBcs
{
public:
    template<class Type, class MeshType>
    static bool match(const boundaryCondition<Type,MeshType>& bc)
    {
        return !bc.eliminated();
    }
};

}

}

}

#endif
