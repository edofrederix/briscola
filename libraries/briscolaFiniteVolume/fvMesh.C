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
    briscola::mesh(dict),
    regIOobject(dict),
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
    return staggeredMetrics_();
}

}

}

}
