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
    immersedBoundaryMethod<Type,MeshType>(dict,fvMsh,false)
{}

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
        if (this->mask_(l,d,i,j,k) == 1)
        {
            // Set sources to 0 in IB
            ls.b()(l,d,i,j,k) = Zero;
        }
    }

    forAllLevels(ls.A(),l,d,i,j,k)
    {
        if (this->mask_(l,d,i,j,k) == 1)
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
