#include "penalizationDirichletImmersedBoundaryCondition.H"
#include "immersedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
penalizationDirichletImmersedBoundaryCondition<Type,MeshType>::penalizationDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>(mshField,ib)
{}

// Destructor

template<class Type, class MeshType>
penalizationDirichletImmersedBoundaryCondition<Type,MeshType>::~penalizationDirichletImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void penalizationDirichletImmersedBoundaryCondition<Type,MeshType>::correctLinearSystem
(
    linearSystem<stencil,Type,MeshType>& ls
)
{
    forAllCells(ls.b(),l,d,i,j,k)
    {
        if (this->IB_.mask()(l,d,i,j,k))
        {
            // Set sources to 0 in IB
            ls.b()(l,d,i,j,k) = Zero;
        }
    }

    forAllCells(ls.A(),l,d,i,j,k)
    {
        if (this->IB_.mask()(l,d,i,j,k))
        {
            // Set coefficients to 0 in IB
            ls.A()(l,d,i,j,k) = Zero;

            // Set central coefficients to 1 in IB
            ls.A()(l,d,i,j,k).center() = 1.0;
        }
    }
}

}

}

}
