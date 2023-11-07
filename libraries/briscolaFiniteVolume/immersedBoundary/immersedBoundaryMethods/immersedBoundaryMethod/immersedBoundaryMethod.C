#include "immersedBoundaryMethod.H"

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
immersedBoundaryMethod<Type,MeshType>::immersedBoundaryMethod
(
    dictionary& dict,
    const fvMesh& fvMsh,
    bool jac
)
:
    fvMsh_(fvMsh),
    massSourceActive_
    (
        dict.lookupOrDefault<bool>("massSourceActive", false)
    ),
    JacobiGhostMethod_(jac),
    omega_
    (
        dict.lookupOrDefault<scalar>("omega", 0.8)
    ),
    ghostMask_
    (
        "ghostMask",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true,
        false,
        true
    )
{
    int nEntries = dict.size();

    // Add shapes to IB according to dictionary entries
    for (int e = 0; e < nEntries; e++)
    {
        word entry = dict.toc()[e];

        if (dict.isDict(entry))
        {
            dictionary entryDict = dict.subDict(entry);

            shapes_.append
            (
                shape::New
                (
                    entryDict,
                    bool(entryDict.lookupOrDefault("inverted", false))
                )
            );
        }
    }

    // Cell centers
    const meshField<vector,MeshType>& CC =
        fvMsh_.metrics<MeshType>().cellCenters();

    // Set IB mask fields
    forAllLevels(ghostMask_,l,d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        ghostMask_(l,d,i,j,k) = 0;

        if (this->isInside(CC(l,d,i,j,k)))
        {
            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];

                if
                (
                    !this->isInside(CC[l][d](ijk+fo))
                )
                {
                    ghostMask_(l,d,i,j,k) = 1;
                }
            }
        }
    }
}

template<class Type, class MeshType>
immersedBoundaryMethod<Type,MeshType>::~immersedBoundaryMethod()
{}

template<class Type, class MeshType>
autoPtr<immersedBoundaryMethod<Type,MeshType>> immersedBoundaryMethod<Type,MeshType>::New
(
    dictionary& dict,
    const fvMesh& fvMsh
)
{
    const word IBMType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(IBMType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown IBM "
            << IBMType << nl << nl
            << "Valid IBM types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<immersedBoundaryMethod<Type,MeshType>>(cstrIter()(dict, fvMsh));
}

template<class Type, class MeshType>
bool immersedBoundaryMethod<Type,MeshType>::isInside(vector xyz)
{
    // Check if xyz is inside any of the IB shapes
    for (int s = 0; s < shapes_.size(); s++)
    {
        if(shapes_[s].isInside(xyz))
        {
            return true;
        }
    }

    return false;
}

template<class Type, class MeshType>
scalar immersedBoundaryMethod<Type,MeshType>::wallDistance(vector c, vector nb)
{
    if (this->isInside(c))
    {
        FatalError
            << "Central point should be fluid node."
            << endl;
        FatalError.exit();
    }

    if (!this->isInside(nb))
    {
        FatalError
            << "Neighbor point should be inside the immersed boundary."
            << endl;
        FatalError.exit();
    }

    scalar dist = -1;

    for (int s = 0; s < shapes_.size(); s++)
    {
        if (mag(dist+1) < 0.01)
        {
            dist = shapes_[s].wallDistance(c, nb);
        }

        if
        (
               (shapes_[s].wallDistance(c, nb) >= 0)
            && (shapes_[s].wallDistance(c, nb) <= dist)
        )
        {
            dist = shapes_[s].wallDistance(c, nb);
        }
    }

    if
    (
           (dist < 0)
        || (dist > mag(c-nb))
    )
    {
        FatalError
            << "No immersed boundary intersection found."
            << " Distance = " << dist << ", mag(c-nb) = " << mag(c-nb)
            << endl;
        FatalError.exit();
    }

    return dist;
}

template<class Type, class MeshType>
scalar immersedBoundaryMethod<Type,MeshType>::wallNormalDistance
(
    vector gc
)
{
    if (!this->isInside(gc))
    {
        FatalError
            << "Ghost cell should be inside the immersed boundary."
            << endl;
        FatalError.exit();
    }

    scalar dist = -1;

    for (int s = 0; s < shapes_.size(); s++)
    {
        if (mag(dist+1) < 0.01)
        {
            dist = shapes_[s].wallNormalDistance(gc);
        }

        if
        (
               (shapes_[s].wallNormalDistance(gc) >= 0)
            && (shapes_[s].wallNormalDistance(gc) <= dist)
        )
        {
            dist = shapes_[s].wallNormalDistance(gc);
        }
    }

    if (dist < 0)
    {
        FatalError
            << "No immersed boundary intersection found."
            << " Distance = " << dist << ", for ghost cell = "
            << gc
            << endl;
        FatalError.exit();
    }

    return dist;
}

template<class Type, class MeshType>
vector immersedBoundaryMethod<Type,MeshType>::mirrorPoint
(
    vector gc
)
{
    if (!this->isInside(gc))
    {
        FatalError
            << "Ghost cell should be inside the immersed boundary."
            << endl;
        FatalError.exit();
    }

    vector mirror = gc;
    scalar dist = -1;

    for (int s = 0; s < shapes_.size(); s++)
    {
        if (mag(dist+1) < 0.01)
        {
            dist = shapes_[s].wallNormalDistance(gc);
            mirror = shapes_[s].mirrorPoint(gc);
        }

        if
        (
               (shapes_[s].wallNormalDistance(gc) >= 0)
            && (shapes_[s].wallNormalDistance(gc) <= dist)
        )
        {
            dist = shapes_[s].wallNormalDistance(gc);
            mirror = shapes_[s].mirrorPoint(gc);
        }
    }

    if (dist < 0)
    {
        FatalError
            << "No immersed boundary intersection found."
            << " Distance = " << dist << ", for ghost cell = "
            << gc
            << endl;
        FatalError.exit();
    }

    return mirror;
}

template<class Type, class MeshType>
void immersedBoundaryMethod<Type,MeshType>::correctLinearSystem
(
    linearSystem<stencil,Type,MeshType>&
)
{
    // Do nothing by default
}

template<class Type, class MeshType>
void immersedBoundaryMethod<Type,MeshType>::correctJacobiPoints
(
    meshLevel<Type,MeshType>&
)
{
    // Do nothing by default
}

template<>
tmp<colocatedScalarField> immersedBoundaryMethod<scalar,staggered>::IBMSource
(
    const staggeredScalarField& field
)
{
    tmp<colocatedScalarField> tSource
    (
        new colocatedScalarField
        (
            "IBMSource",
            field.fvMsh()
        )
    );

    colocatedScalarField& source = tSource.ref();

    source = Zero;

    if (JacobiGhostMethod_)
    {
        const colocatedFaceScalarField& fa =
            fvMsh_.metrics<colocated>().faceAreas();

        const colocatedScalarField& cv =
            fvMsh_.metrics<colocated>().cellVolumes();

        forAllCells(source[0][0],i,j,k)
        {
            source(0,0,i,j,k) -=
                ghostMask_(0,0,i,j,k) * field(0,0,i,j,k) * fa(0,0,i,j,k).left();

            source(0,0,i,j,k) +=
                ghostMask_(0,0,i+1,j,k) * field(0,0,i+1,j,k) * fa(0,0,i,j,k).right();

            source(0,0,i,j,k) -=
                ghostMask_(0,1,i,j,k) * field(0,1,i,j,k) * fa(0,0,i,j,k).bottom();

            source(0,0,i,j,k) +=
                ghostMask_(0,1,i,j+1,k) * field(0,1,i,j+1,k) * fa(0,0,i,j,k).top();

            source(0,0,i,j,k) -=
                ghostMask_(0,2,i,j,k) * field(0,2,i,j,k) * fa(0,0,i,j,k).aft();

            source(0,0,i,j,k) +=
                ghostMask_(0,2,i,j,k+1) * field(0,2,i,j,k+1) * fa(0,0,i,j,k).fore();

            source(0,0,i,j,k) /= cv(0,0,i,j,k);
        }
    }

    return tSource;
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam
