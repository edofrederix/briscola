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

template<class MeshType>
void fvMesh::setInternalCells()
{
    forAll(*this, l)
    {
        const partLevel& level = this->operator[](l);

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            faceLabel& I = this->I<MeshType>(l,d);

            const labelVector& padding = MeshType::padding[d];
            const faceLabel slave = mshPtr_->facePatchSlave();

            I = faceLabel(zeroXYZ, level.N()+padding);

            for (int i = 0; i < 6; i++)
                I[i] += (padding[i/2] && slave[i]) ? 1 - 2*(i%2) : 0;
        }
    }
}

template<>
faceLabel& fvMesh::I<colocated>(const label l, const label d)
{
    return Ic_[l][d];
}

template<>
faceLabel& fvMesh::I<staggered>(const label l, const label d)
{
    return Is_[l][d];
}

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
    Ic_(mshPtr_->size(), List<faceLabel>(colocated::numberOfDirections)),
    Is_(mshPtr_->size(), List<faceLabel>(staggered::numberOfDirections)),
    colocatedMetrics_(),
    staggeredMetrics_()
{
    setInternalCells<colocated>();
    setInternalCells<staggered>();

    colocatedMetrics_ = new fvMeshMetrics<colocated>(*this);

    // Only generate staggered metrics when the brick topology is structured

    if (mshPtr_->structured())
        staggeredMetrics_ = new fvMeshMetrics<staggered>(*this);
}

fvMesh::fvMesh(const fvMesh& fvMsh)
:
    regIOobject(fvMsh, true),
    mshPtr_(fvMsh.mshPtr_, false),
    schemeDict_(fvMsh.schemeDict_),
    solverDict_(fvMsh.solverDict_),
    Ic_(fvMsh.Ic_),
    Is_(fvMsh.Is_),
    colocatedMetrics_(),
    staggeredMetrics_()
{
    colocatedMetrics_ = new fvMeshMetrics<colocated>(*this);

    // Only generate staggered metrics when the brick topology is structured

    if (mshPtr_->structured())
        staggeredMetrics_ = new fvMeshMetrics<staggered>(*this);
}

fvMesh::~fvMesh()
{}

template<>
const faceLabel& fvMesh::I<colocated>(const label l, const label d) const
{
    return Ic_[l][d];
}

template<>
const faceLabel& fvMesh::I<staggered>(const label l, const label d) const
{
    return Is_[l][d];
}

template<>
const fvMeshMetrics<colocated>& fvMesh::metrics<colocated>() const
{
    return colocatedMetrics_();
}

template<>
const fvMeshMetrics<staggered>& fvMesh::metrics<staggered>() const
{
    #ifdef FULLDEBUG
    if (!mshPtr_->structured())
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
