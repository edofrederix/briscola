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
penalizationDirichletImmersedBoundaryCondition<Type,MeshType>::
penalizationDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>
    (
        mshField,
        ib,
        ib.mask()
    ),
    boundaryValues_(this->dict().lookup("values"))
{}

// Destructor

template<class Type, class MeshType>
penalizationDirichletImmersedBoundaryCondition<Type,MeshType>::
~penalizationDirichletImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void penalizationDirichletImmersedBoundaryCondition<Type,MeshType>
::correctJacobiPoints
(
    meshLevel<Type,MeshType>& x
) const
{
    scalar omega = this->omega_;

    label l = x.levelNum();

    if (l == 0)
    {
        forAllCells(x,d,i,j,k)
        {
            if (this->forcingPoints_(l,d,i,j,k))
            {
                x(d,i,j,k) = (1.0 - omega) * x(d,i,j,k)
                    + omega * boundaryValues_[d];
            }
        }
    }
}

}

}

}
