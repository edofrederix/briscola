#include "immersedBoundaryCondition.H"
#include "fvMesh.H"
#include "meshField.H"
#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"
#include "linearSystem.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

using Foam::max;
using Foam::min;
using Foam::sqr;
using Foam::mag;

// Constructor

template<class Type, class MeshType>
immersedBoundaryCondition<Type,MeshType>::immersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    fvMsh_(mshField.fvMsh()),
    IB_(ib),
    dict_
    (
        mshField.found("boundaryConditions")
     && mshField.subDict("boundaryConditions").found("ImmersedBoundary")
      ? mshField.subDict("boundaryConditions").subDict("ImmersedBoundary")
      : dictionary::null
    )
{}

template<class Type, class MeshType>
immersedBoundaryCondition<Type,MeshType>::~immersedBoundaryCondition()
{}

template<class Type, class MeshType>
autoPtr<immersedBoundaryCondition<Type,MeshType>>
immersedBoundaryCondition<Type,MeshType>::New
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
{
    dictionary IBCDict
    (
        mshField.subDict("boundaryConditions").subDict("ImmersedBoundary")
    );

    const word IBCType(IBCDict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(IBCType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown immersed boundary condition "
            << IBCType << nl << nl
            << "Valid immersed boundary condition types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<immersedBoundaryCondition<Type,MeshType>>(cstrIter()(mshField, ib));
}

template<class Type, class MeshType>
void immersedBoundaryCondition<Type,MeshType>::correctLinearSystem
(
    linearSystem<stencil,Type,MeshType>&
)
{
    // Do nothing by default
}

template<class Type, class MeshType>
void immersedBoundaryCondition<Type,MeshType>::correctJacobiPoints
(
    meshLevel<Type,MeshType>&
)
{
    // Do nothing by default
}

// template<>
// tmp<colocatedScalarField> immersedBoundaryCondition<scalar,staggered>::IBMSource
// (
//     const staggeredScalarField& field
// )
// {
//     tmp<colocatedScalarField> tSource
//     (
//         new colocatedScalarField
//         (
//             "IBMSource",
//             field.fvMsh()
//         )
//     );

//     colocatedScalarField& source = tSource.ref();

//     source = Zero;

//     if (JacobiGhostMethod_)
//     {
//         const colocatedFaceScalarField& fa =
//             fvMsh_.metrics<colocated>().faceAreas();

//         const colocatedScalarField& cv =
//             fvMsh_.metrics<colocated>().cellVolumes();

//         forAllCells(source[0][0],i,j,k)
//         {
//             source(0,0,i,j,k) -=
//                 ghostMask_(0,0,i,j,k) * field(0,0,i,j,k) * fa(0,0,i,j,k).left();

//             source(0,0,i,j,k) +=
//                 ghostMask_(0,0,i+1,j,k) * field(0,0,i+1,j,k) * fa(0,0,i,j,k).right();

//             source(0,0,i,j,k) -=
//                 ghostMask_(0,1,i,j,k) * field(0,1,i,j,k) * fa(0,0,i,j,k).bottom();

//             source(0,0,i,j,k) +=
//                 ghostMask_(0,1,i,j+1,k) * field(0,1,i,j+1,k) * fa(0,0,i,j,k).top();

//             source(0,0,i,j,k) -=
//                 ghostMask_(0,2,i,j,k) * field(0,2,i,j,k) * fa(0,0,i,j,k).aft();

//             source(0,0,i,j,k) +=
//                 ghostMask_(0,2,i,j,k+1) * field(0,2,i,j,k+1) * fa(0,0,i,j,k).fore();

//             source(0,0,i,j,k) /= cv(0,0,i,j,k);
//         }
//     }

//     return tSource;
// }

}

}

}
