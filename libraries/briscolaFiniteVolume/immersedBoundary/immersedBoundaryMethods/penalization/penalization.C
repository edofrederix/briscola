#include "penalization.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
penalization<Type,MeshType>::penalization
(
    dictionary& dict,
    const fvMesh& fvMsh
)
:
    immersedBoundaryMethod<Type,MeshType>(dict,fvMsh,false),
    mask_
    (
        "mask",
        fvMsh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true,
        false,
        true
    )
{
    // Cell centers
    const meshField<vector,MeshType>& CC =
        fvMsh.metrics<MeshType>().cellCenters();

    // Set IB mask fields
    forAllLevels(mask_,l,d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        mask_(l,d,i,j,k) = 0;

        if (this->isInside(CC(l,d,i,j,k)))
        {
            mask_(l,d,i,j,k) = 1;
        }
    }
}

// Destructor

template<class Type, class MeshType>
penalization<Type,MeshType>::~penalization()
{}

template<class Type, class MeshType>
void penalization<Type,MeshType>::correctLinearSystem
(
    linearSystem<stencil,Type,MeshType>& ls
)
{
    forAllLevels(ls.b(),l,d,i,j,k)
    {
        if (mask_(l,d,i,j,k) == 1)
        {
            // Set sources to 0 in IB
            ls.b()(l,d,i,j,k) = Zero;
        }
    }

    forAllLevels(ls.A(),l,d,i,j,k)
    {
        if (mask_(l,d,i,j,k) == 1)
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
