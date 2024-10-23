#include "fvMesh.H"

#include "colocatedFields.H"
#include "staggeredFields.H"
#include "immersedBoundary.H"

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
        const part& p = this->operator[](l);

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            faceLabel& I = this->I<MeshType>(l,d);

            const labelVector& padding = MeshType::padding[d];
            const faceLabel slave = mshPtr_->facePatchSlave();

            I = faceLabel(zeroXYZ, p.N()+padding);

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
    {
        staggeredMetrics_ = new fvMeshMetrics<staggered>(*this);
    }
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
    colocatedMetrics_->setIBs();

    // Only generate staggered metrics when the brick topology is structured

    if (mshPtr_->structured())
    {
        staggeredMetrics_ = new fvMeshMetrics<staggered>(*this);
        staggeredMetrics_->setIBs();
    }
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

template<>
labelVector fvMesh::findCell<colocated>
(
    const vector& p,
    const label l,
    const label
) const
{
    return mshPtr_->findCell(p,l);
}

template<>
labelVector fvMesh::findCell<staggered>
(
    const vector& p,
    const label l,
    const label d
) const
{
    labelVector colo = mshPtr_->findCell(p,l);

    const meshDirection<vertexVector,staggered>& v =
        this->metrics<staggered>().vertexCenters()[l][d];

    labelVector vU(v.I().upper());

    if (colo != -unitXYZ)
    {
        // Search staggered cells that may overlap with the colocated one

        labelVector L(colo);
        labelVector U(colo + staggered::padding[d] + unitXYZ);

        if (this->rectilinear() != unitXYZ)
        {
            // Add surrounding cells in non-shifted direction

            for (int dir = 0; dir < 3; dir++)
            if (dir != d)
            {
                L -= staggered::padding[dir];
                U += staggered::padding[dir];
            }

            // Trim

            for (int dir = 0; dir < 3; dir++)
            if (dir != d)
            {
                L[dir] = Foam::max(L[dir], 0);
                U[dir] = Foam::min(U[dir], vU[dir]);
            }
        }

        for (int i = L.x(); i < U.x(); i++)
            for (int j = L.y(); j < U.y(); j++)
                for (int k = L.z(); k < U.z(); k++)
                    if (interpolationWeights(p,v(i,j,k),true) != -vector::one)
                        return labelVector(i,j,k);

        // When the point is near a shifted boundary, it may be found on the
        // colocated mesh but not on the staggered mesh.

        return -unitXYZ;
    }
    else
    {
        // Near shifted boundaries, the point may be outside the colocated mesh
        // but still on the staggered mesh. This happens in two situations:
        //
        // 1) The point is outside the physical domain -> return -unitXYZ
        // 2) The point is on a parallel/periodic master/slave slab of cells.
        //    Data on such cells is redundantly stored across two processors, so
        //    we let the processor which has the point on its colocated mesh
        //    handle the situation -> return -unitXYZ

        return -unitXYZ;
    }
}

template<class MeshType>
List<labelVector> fvMesh::findCells
(
    const vectorList& points,
    const label l,
    const label d
) const
{
    List<labelVector> res(points.size());

    forAll(points, i)
        res[i] = this->findCell<MeshType>(points[i], l, d);

    return res;
}

// Instantiate

template List<labelVector>
fvMesh::findCells<colocated>(const vectorList&, const label, const label) const;
template List<labelVector>
fvMesh::findCells<staggered>(const vectorList&, const label, const label) const;

template const PtrList<immersedBoundary<staggered>>&
fvMesh::IBs<staggered>() const;
template const PtrList<immersedBoundary<colocated>>&
fvMesh::IBs<colocated>() const;

}

}

}
