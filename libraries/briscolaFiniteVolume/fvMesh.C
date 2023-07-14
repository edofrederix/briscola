#include "fvMesh.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(fvMesh, 0);

fvMesh::fvMesh(const IOdictionary& dict, const Time& time)
:
    regIOobject(dict, true),
    mshPtr_(mesh::New(dict)),
    schemeDict_
    (
        IOobject
        (
            "briscolaSchemeDict",
            time.system(),
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    solverDict_
    (
        IOobject
        (
            "briscolaSolverDict",
            time.system(),
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    colocatedMetrics_(),
    staggeredMetrics_()
{
    colocatedMetrics_ = new fvMeshMetrics<colocated>(*this);

    // Only generate staggered metrics when the brick topology is structured

    if (mshPtr_->topology().structured())
        staggeredMetrics_ = new fvMeshMetrics<staggered>(*this);
}

fvMesh::fvMesh(const fvMesh& fvMsh)
:
    regIOobject(fvMsh, true),
    mshPtr_(fvMsh.mshPtr_, false),
    schemeDict_(fvMsh.schemeDict_),
    solverDict_(fvMsh.solverDict_),
    colocatedMetrics_(),
    staggeredMetrics_()
{
    colocatedMetrics_ = new fvMeshMetrics<colocated>(*this);

    // Only generate staggered metrics when the brick topology is structured

    if (mshPtr_->topology().structured())
        staggeredMetrics_ = new fvMeshMetrics<staggered>(*this);
}

fvMesh::~fvMesh()
{}

template<>
const fvMeshMetrics<colocated>& fvMesh::metrics<colocated>() const
{
    return colocatedMetrics_();
}

template<>
const fvMeshMetrics<staggered>& fvMesh::metrics<staggered>() const
{
    #if DEBUG
    if (!topology().structured())
    {
        FatalErrorInFunction
            << "Staggered metrics are not generated on unstructured meshes."
            << endl << abort(FatalError);
    }
    #endif

    return staggeredMetrics_();
}

}

}

}
