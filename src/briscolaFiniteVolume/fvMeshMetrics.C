#include "fvMeshMetrics.H"

#include "colocatedFields.H"
#include "staggeredFields.H"
#include "faceFields.H"

#include "fvMesh.H"
#include "immersedBoundary.H"

#include "geometricObjects.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
void fvMeshMetrics<MeshType>::setCellCenters()
{
    const meshField<vertexVector,MeshType>& vc = vertexCenters();

    cellCentersPtr_ =
        new meshField<vector,MeshType>
        (
            word(MeshType::typeName) + "CellCenters",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    meshField<vector,MeshType>& cc = *cellCentersPtr_;

    cc = Zero;

    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                for (int vi = 0; vi < 8; vi++)
                    cc(l,d,ijk) = hexa(vc(l,d,ijk)).center();
            }
        }
    }

    cc.correctAggData();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setCellVolumes()
{
    const meshField<vertexVector,MeshType>& vc = vertexCenters();

    cellVolumesPtr_ =
        new meshField<scalar,MeshType>
        (
            word(MeshType::typeName) + "CellVolumes",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    meshField<scalar,MeshType>& cv = *cellVolumesPtr_;

    cv = Zero;

    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            // Cell volume computed from tet decomposition, needed for
            // consistency with geometric VoF calculations.

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                cv(l,d,ijk) = hexa(vc(l,d,ijk)).volume();
            }
        }
    }

    // Small number to avoid division by zero in unset/invalid ghost cells
    cv = max(cv, 1e-16);
    cv.correctAggData();

    // Compute at store the inverse cell volumes

    inverseCellVolumesPtr_ =
        new meshField<scalar,MeshType>
        (
            word(MeshType::typeName) + "InverseCellVolumes",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    meshField<scalar,MeshType>& icv = *inverseCellVolumesPtr_;

    icv = 1.0/cv;
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setFaceCenters()
{
    const meshField<faceVector,MeshType> fcAoS(aos().faceCenters());

    faceCentersPtr_ =
        new faceField<vector,MeshType>
        (
            word(MeshType::typeName) + "FaceCenters",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    faceField<vector,MeshType>& fc = *faceCentersPtr_;

    fc = Zero;

    // Use a cell loop that also traverses ghost cells. We cannot use a boundary
    // correction on SoA storage because face directions may not be aligned
    // across processor boundaries.

    for (int fd = 0; fd < 3; fd++)
    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                fc[fd](l,d,ijk) = fcAoS(l,d,ijk)[fd*2];
            }
        }
    }

    for (int fd = 0; fd < 3; fd++)
        fc[fd].correctAggData();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setFaceNormals()
{
    const meshField<faceVector,MeshType> fnAoS(aos().faceNormals());

    faceNormalsPtr_ =
        new faceField<vector,MeshType>
        (
            word(MeshType::typeName) + "FaceNormals",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    faceField<vector,MeshType>& fn = *faceNormalsPtr_;

    fn = Zero;

    // Use a cell loop that also traverses ghost cells. We cannot use a boundary
    // correction on SoA storage because face directions may not be aligned
    // across processor boundaries.

    for (int fd = 0; fd < 3; fd++)
    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                fn[fd](l,d,ijk) = fnAoS(l,d,ijk)[fd*2];
            }
        }
    }

    for (int fd = 0; fd < 3; fd++)
        fn[fd].correctAggData();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setFaceAreas()
{
    const meshField<faceScalar,MeshType> faAoS(aos().faceAreas());

    faceAreasPtr_ =
        new faceField<scalar,MeshType>
        (
            word(MeshType::typeName) + "FaceAreas",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    faceField<scalar,MeshType>& fa = *faceAreasPtr_;

    fa = Zero;

    // Use a cell loop that also traverses ghost cells. We cannot use a boundary
    // correction on SoA storage because face directions may not be aligned
    // across processor boundaries.

    for (int fd = 0; fd < 3; fd++)
    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                fa[fd](l,d,ijk) = faAoS(l,d,ijk)[fd*2];
            }
        }
    }

    for (int fd = 0; fd < 3; fd++)
        fa[fd].correctAggData();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setFaceAreaNormals()
{
    const meshField<faceVector,MeshType> fanAoS(aos().faceAreaNormals());

    faceAreaNormalsPtr_ =
        new faceField<vector,MeshType>
        (
            word(MeshType::typeName) + "FaceAreaNormals",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    faceField<vector,MeshType>& fan = *faceAreaNormalsPtr_;

    fan = Zero;

    // Use a cell loop that also traverses ghost cells. We cannot use a boundary
    // correction on SoA storage because face directions may not be aligned
    // across processor boundaries.

    for (int fd = 0; fd < 3; fd++)
    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                fan[fd](l,d,ijk) = fanAoS(l,d,ijk)[fd*2];
            }
        }
    }

    for (int fd = 0; fd < 3; fd++)
        fan[fd].correctAggData();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setFaceDeltas()
{
    const meshField<faceScalar,MeshType> deltaAoS(aos().faceDeltas());

    faceDeltasPtr_ =
        new faceField<scalar,MeshType>
        (
            word(MeshType::typeName) + "FaceDeltas",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    faceField<scalar,MeshType>& delta = *faceDeltasPtr_;

    delta = Zero;

    // Use a cell loop that also traverses ghost cells. We cannot use a boundary
    // correction on SoA storage because face directions may not be aligned
    // across processor boundaries.

    for (int fd = 0; fd < 3; fd++)
    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                delta[fd](l,d,ijk) = deltaAoS(l,d,ijk)[fd*2];
            }
        }
    }

    for (int fd = 0; fd < 3; fd++)
        delta[fd].correctAggData();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setFaceWeights()
{
    const meshField<faceScalar,MeshType> fwcAoS(aos().faceWeightsCenter());
    const meshField<faceScalar,MeshType> fwnAoS(aos().faceWeightsNeighbor());

    faceWeightsCenterPtr_ =
        new faceField<scalar,MeshType>
        (
            word(MeshType::typeName) + "FaceWeightsCenter",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    faceWeightsNeighborPtr_ =
        new faceField<scalar,MeshType>
        (
            word(MeshType::typeName) + "FaceWeightsNeighbor",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    faceField<scalar,MeshType>& fwc = *faceWeightsCenterPtr_;
    faceField<scalar,MeshType>& fwn = *faceWeightsNeighborPtr_;

    // Initialize weights as 50/50

    fwc = 0.5;
    fwn = 0.5;

    // Use a cell loop that also traverses ghost cells. We cannot use a boundary
    // correction on SoA storage because face directions may not be aligned
    // across processor boundaries.

    for (int fd = 0; fd < 3; fd++)
    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                fwc[fd](l,d,ijk) = fwcAoS(l,d,ijk)[fd*2];
                fwn[fd](l,d,ijk) = fwnAoS(l,d,ijk)[fd*2];
            }
        }
    }

    for (int fd = 0; fd < 3; fd++)
    {
        fwc[fd].correctAggData();
        fwn[fd].correctAggData();
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setVertexCenters()
{
    vertexCentersPtr_ =
        new meshField<vertexVector,MeshType>
        (
            word(MeshType::typeName) + "VertexCenters",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    meshField<vertexVector,MeshType>& vc = *vertexCentersPtr_;

    vc = Zero;

    forAll(fvMsh_, l)
    {
        const levelPoints& points = fvMsh_[l].points();

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);
            const vector shift = MeshType::shift[d];

            // Create vertices from a levelPoints object, which consistently
            // generates ghost points

            levelPoints p(fvMsh_[l]);

            p.clear();
            p.setSizeFromCells(N);

            forAllBlock(p, i, j, k)
                p(i,j,k) = points.interp(vector(i,j,k) + shift);

            p.calcGhostPoints();
            p.clean();

            // Nothing to do on empty levels

            if (fvMsh_[l].empty())
                continue;

            // For each cell, copy points into vertex vectors

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                vc(l,d,ijk) =
                    vertexVector
                    (
                        p(ijk),
                        p(ijk + unitX),
                        p(ijk + unitY),
                        p(ijk + unitXY),
                        p(ijk + unitZ),
                        p(ijk + unitXZ),
                        p(ijk + unitYZ),
                        p(ijk + unitXYZ)
                    );

                // Create hexahedral from vertex vector and check its validity.
                // On unstructured meshes only internal cells and face ghosts
                // need to be valid. Subsequent functions do not need to check
                // anymore.

                const hexa h(vc(l,d,ijk));

                if
                (
                    fvMsh_.structured()
                 || (
                      - cmptSum(briscola::cmptMin(ijk,zeroXYZ))
                      + cmptSum(briscola::cmptMax(ijk-N+unitXYZ,zeroXYZ))
                     <= 1
                    )
                )
                    if (!h.valid())
                        FatalErrorInFunction
                            << "Invalid hexahedral created from vertices "
                            << vc(l,d,ijk) << nl << "at"
                            << " level " << l
                            << " direction " << d
                            << " index " << ijk << endl
                            << abort(FatalError);
            }
        }
    }

    vc.correctAggData();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setGlobalCellNumbers()
{
    globalCellNumbersPtr_ =
        new meshField<label,MeshType>
        (
            word(MeshType::typeName) + "GlobalCellNumbers",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    meshField<label,MeshType>& gcn = *globalCellNumbersPtr_;
    gcn = -1;

    globalCellCounts_.resize(gcn.size());

    // Cells are numbered lexicographically on each processor and with
    // increasing processor rank

    forAll(gcn, l)
    {
        globalCellCounts_[l].resize(MeshType::numberOfDirections);

        forAll(gcn[l], d)
        {
            labelList sizes(Pstream::nProcs());

            sizes[Pstream::myProcNo()] =
                cmptProduct(fvMsh_.N<MeshType>(l,d));

            Pstream::gatherList(sizes, Pstream::msgType(), fvMsh_[l].comms());
            Pstream::scatterList(sizes, Pstream::msgType(), fvMsh_[l].comms());

            label start = 0;

            for (int i = 0; i < Pstream::myProcNo(); i++)
                start += sizes[i];

            int c = 0;
            forAllCells(gcn[l][d], i, j, k)
                gcn[l][d](i,j,k) = start + c++;

            globalCellCounts_[l][d] = sum(sizes);
        }
    }

    gcn.correctCommsBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setImmersedBoundaries()
{
    if (fvMsh_.msh().dict().found("immersedBoundaries"))
    {
        const dictionary& ibmDict =
            fvMsh_.msh().dict().subDict("immersedBoundaries");

        label size = ibmDict.size();

        for (int i = 0; i < size; i++)
        {
            word IBname = ibmDict.toc()[i];

            immersedBoundaries_.append
            (
                new immersedBoundary<MeshType>
                (
                    *this,
                    ibmDict.subDict(IBname)
                )
            );
        }
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setImmersedBoundaryMask()
{
    maskPtr_ =
        new meshField<scalar,MeshType>
        (
            word(MeshType::typeName) + "Mask",
            fvMsh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true,
            true
        );

    meshField<scalar,MeshType>& mask = *maskPtr_;
    mask = Zero;

    if (immersedBoundaries_.size())
    {
        forAll(immersedBoundaries_, i)
            mask += immersedBoundaries_[i].mask();

        mask = min(mask, 1.0);
    }
}

template<class MeshType>
fvMeshMetrics<MeshType>::fvMeshMetrics(const fvMesh& fvMsh)
:
    fvMsh_(fvMsh),
    cellCentersPtr_(nullptr),
    cellVolumesPtr_(nullptr),
    inverseCellVolumesPtr_(nullptr),
    faceCentersPtr_(nullptr),
    faceNormalsPtr_(nullptr),
    faceAreasPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    faceDeltasPtr_(nullptr),
    faceWeightsCenterPtr_(nullptr),
    faceWeightsNeighborPtr_(nullptr),
    vertexCentersPtr_(nullptr),
    globalCellNumbersPtr_(nullptr),
    maskPtr_(nullptr),
    aosPtr_(new AoS(*this))
{
    // These need to be set a priori
    setGlobalCellNumbers();
    setImmersedBoundaries();
}

template<class MeshType>
fvMeshMetrics<MeshType>::~fvMeshMetrics()
{
    // Cells

    if (cellCentersPtr_)
        delete cellCentersPtr_;

    if (cellVolumesPtr_)
        delete cellVolumesPtr_;

    if (inverseCellVolumesPtr_)
        delete inverseCellVolumesPtr_;

    // Faces

    if (faceCentersPtr_)
        delete faceCentersPtr_;

    if (faceNormalsPtr_)
        delete faceNormalsPtr_;

    if (faceAreasPtr_)
        delete faceAreasPtr_;

    if (faceAreaNormalsPtr_)
        delete faceAreaNormalsPtr_;

    if (faceDeltasPtr_)
        delete faceDeltasPtr_;

    if (faceWeightsCenterPtr_)
        delete faceWeightsCenterPtr_;

    if (faceWeightsNeighborPtr_)
        delete faceWeightsNeighborPtr_;

    // Vertices

    if (vertexCentersPtr_)
        delete vertexCentersPtr_;

    // Cell numbers

    if (globalCellNumbersPtr_)
        delete globalCellNumbersPtr_;

    // Mask

    if (maskPtr_)
        delete maskPtr_;
}

// On-demand Array of Structures (AoS) metrics class

template<class MeshType>
tmp<meshField<faceVector,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceCenters() const
{
    const meshField<vertexVector,MeshType>& vc = metrics_.vertexCenters();

    tmp<meshField<faceVector,MeshType>> tFc =
        meshField<faceVector,MeshType>::New("faceCenters", metrics_.fvMsh_);

    meshField<faceVector,MeshType>& fc = tFc.ref();
    fc.makeDeep();

    fc = Zero;

    forAll(metrics_.fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = metrics_.fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                const hexa h(vc(l,d,ijk));

                for (int fi = 0; fi < 6; fi++)
                    fc(l,d,ijk)[fi] = h.faceCenter(fi);
            }
        }
    }

    fc.correctAggData();

    return tFc;
}

template<class MeshType>
tmp<meshField<faceVector,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceNormals() const
{
    const meshField<vertexVector,MeshType>& vc = metrics_.vertexCenters();

    tmp<meshField<faceVector,MeshType>> tFn =
        meshField<faceVector,MeshType>::New("faceNormals", metrics_.fvMsh_);

    meshField<faceVector,MeshType>& fn = tFn.ref();
    fn.makeDeep();

    fn = Zero;

    forAll(metrics_.fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = metrics_.fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                for (int fi = 0; fi < 6; fi++)
                {
                    fn(l,d,ijk)[fi] =
                        Foam::normalised(hexa(vc(l,d,ijk)).faceAreaNormal(fi));
                }
            }
        }
    }

    fn.correctAggData();

    return tFn;
}

template<class MeshType>
tmp<meshField<faceScalar,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceAreas() const
{
    const meshField<vertexVector,MeshType>& vc = metrics_.vertexCenters();

    tmp<meshField<faceScalar,MeshType>> tFa =
        meshField<faceScalar,MeshType>::New("faceAreas", metrics_.fvMsh_);

    meshField<faceScalar,MeshType>& fa = tFa.ref();
    fa.makeDeep();

    fa = Zero;

    forAll(metrics_.fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = metrics_.fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                for (int fi = 0; fi < 6; fi++)
                {
                    fa(l,d,ijk)[fi] =
                        Foam::mag(hexa(vc(l,d,ijk)).faceAreaNormal(fi));
                }
            }
        }
    }

    fa.correctAggData();

    return tFa;
}

template<class MeshType>
tmp<meshField<faceVector,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceAreaNormals() const
{
    const meshField<vertexVector,MeshType>& vc = metrics_.vertexCenters();

    tmp<meshField<faceVector,MeshType>> tFan =
        meshField<faceVector,MeshType>::New("faceAreaNormals", metrics_.fvMsh_);

    meshField<faceVector,MeshType>& fan = tFan.ref();
    fan.makeDeep();

    fan = Zero;

    forAll(metrics_.fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = metrics_.fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                for (int fi = 0; fi < 6; fi++)
                {
                    fan(l,d,ijk)[fi] = hexa(vc(l,d,ijk)).faceAreaNormal(fi);
                }
            }
        }
    }

    fan.correctAggData();

    return tFan;
}

template<class MeshType>
tmp<meshField<faceScalar,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceDeltas() const
{
    const meshField<vector,MeshType>& cc = metrics_.cellCenters();

    tmp<meshField<faceScalar,MeshType>> tDelta =
        meshField<faceScalar,MeshType>::New("faceDeltas", metrics_.fvMsh_);

    meshField<faceScalar,MeshType>& delta = tDelta.ref();
    delta.makeDeep();

    delta = Zero;

    forAll(metrics_.fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = metrics_.fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = 0; ijk.x() < N.x(); ijk.x()++)
            for (ijk.y() = 0; ijk.y() < N.y(); ijk.y()++)
            for (ijk.z() = 0; ijk.z() < N.z(); ijk.z()++)
            {
                for (label f = 0; f < 6; f++)
                {
                    const labelVector nei(neighbor(ijk,f));

                    delta(l,d,ijk)[f] =
                        1.0/Foam::mag(cc(l,d,ijk) - cc(l,d,nei));

                    // Assure that ghost cells also have the relevant face
                    // values set

                    delta(l,d,nei)[f%2 ? f-1 : f+1] = delta(l,d,ijk)[f];
                }
            }
        }
    }

    // Also set remaining ghost cell face values at processor interfaces

    delta.correctCommsBoundaryConditions();
    delta.correctAggData();

    return tDelta;
}

template<class MeshType>
tmp<meshField<faceScalar,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceWeightsCenter() const
{
    const meshField<vector,MeshType>& cc = metrics_.cellCenters();

    const meshField<faceScalar,MeshType> delta(faceDeltas());
    const meshField<faceVector,MeshType> fc(faceCenters());

    tmp<meshField<faceScalar,MeshType>> tFwc =
        meshField<faceScalar,MeshType>::New
        (
            "faceWeightsCenter",
            metrics_.fvMsh_
        );

    meshField<faceScalar,MeshType>& fwc = tFwc.ref();
    fwc.makeDeep();

    fwc = 0.5;

    forAll(metrics_.fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = metrics_.fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = 0; ijk.x() < N.x(); ijk.x()++)
            for (ijk.y() = 0; ijk.y() < N.y(); ijk.y()++)
            for (ijk.z() = 0; ijk.z() < N.z(); ijk.z()++)
            {
                for (label f = 0; f < 6; f++)
                {
                    const labelVector nei(neighbor(ijk,f));

                    // Project the vector that connects the face center and
                    // the cell center onto the cell-to-cell vector

                    fwc(l,d,ijk)[f] =
                        Foam::mag
                        (
                            (cc(l,d,nei) - fc(l,d,ijk)[f])
                          & (cc(l,d,nei) - cc(l,d,ijk))
                        )
                      * Foam::sqr(delta(l,d,ijk)[f]);

                    const scalar fwn =
                        Foam::mag
                        (
                            (cc(l,d,ijk) - fc(l,d,ijk)[f])
                          & (cc(l,d,nei) - cc(l,d,ijk))
                        )
                      * Foam::sqr(delta(l,d,ijk)[f]);

                    // Assure that ghost cells also have the relevant face
                    // values set

                    fwc(l,d,nei)[f%2 ? f-1 : f+1] = fwn;

                    // Check if the weights add up to unity up to some precision

                    const scalar sum =
                        round((fwc(l,d,ijk)[f] + fwn)*1e12)/1e12;

                    if (sum != 1.0)
                        FatalErrorInFunction
                            << "Face weights do not add up to unity at "
                            << "level " << l
                            << ", " << MeshType::typeName << " direction " << d
                            << ", cell " << ijk
                            << ", face " << f << endl << abort(FatalError);
                }
            }
        }
    }

    // Also set remaining ghost cell face values at processor interfaces

    fwc.correctCommsBoundaryConditions();
    fwc.correctAggData();

    return tFwc;
}

template<class MeshType>
tmp<meshField<faceScalar,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceWeightsNeighbor() const
{
    const meshField<vector,MeshType>& cc = metrics_.cellCenters();

    const meshField<faceScalar,MeshType> delta(faceDeltas());
    const meshField<faceVector,MeshType> fc(faceCenters());

    tmp<meshField<faceScalar,MeshType>> tFwn =
        meshField<faceScalar,MeshType>::New
        (
            "faceWeightsNeighbor",
            metrics_.fvMsh_
        );

    meshField<faceScalar,MeshType>& fwn = tFwn.ref();
    fwn.makeDeep();

    fwn = 0.5;

    forAll(metrics_.fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = metrics_.fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = 0; ijk.x() < N.x(); ijk.x()++)
            for (ijk.y() = 0; ijk.y() < N.y(); ijk.y()++)
            for (ijk.z() = 0; ijk.z() < N.z(); ijk.z()++)
            {
                for (label f = 0; f < 6; f++)
                {
                    const labelVector nei(neighbor(ijk,f));

                    // Project the vector that connects the face center and
                    // the cell center onto the cell-to-cell vector

                    const scalar fwc =
                        Foam::mag
                        (
                            (cc(l,d,nei) - fc(l,d,ijk)[f])
                          & (cc(l,d,nei) - cc(l,d,ijk))
                        )
                      * Foam::sqr(delta(l,d,ijk)[f]);

                    fwn(l,d,ijk)[f] =
                        Foam::mag
                        (
                            (cc(l,d,ijk) - fc(l,d,ijk)[f])
                          & (cc(l,d,nei) - cc(l,d,ijk))
                        )
                      * Foam::sqr(delta(l,d,ijk)[f]);

                    // Assure that ghost cells also have the relevant face
                    // values set

                    fwn(l,d,nei)[f%2 ? f-1 : f+1] = fwc;

                    // Check if the weights add up to unity up to some precision

                    const scalar sum =
                        round((fwc + fwn(l,d,ijk)[f])*1e12)/1e12;

                    if (sum != 1.0)
                        FatalErrorInFunction
                            << "Face weights do not add up to unity at "
                            << "level " << l
                            << ", " << MeshType::typeName << " direction " << d
                            << ", cell " << ijk
                            << ", face " << f << endl << abort(FatalError);
                }
            }
        }
    }

    // Also set remaining ghost cell face values at processor interfaces

    fwn.correctCommsBoundaryConditions();
    fwn.correctAggData();

    return tFwn;
}

template<class MeshType>
tmp<meshField<edgeVector,MeshType>>
fvMeshMetrics<MeshType>::AoS::edgeCenters() const
{
    const meshField<vertexVector,MeshType>& vc = metrics_.vertexCenters();

    tmp<meshField<edgeVector,MeshType>> tEc =
        meshField<edgeVector,MeshType>::New("edgeCenters", metrics_.fvMsh_);

    meshField<edgeVector,MeshType>& ec = tEc.ref();
    ec.makeDeep();

    ec = Zero;

    forAll(metrics_.fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = metrics_.fvMsh_.N<MeshType>(l,d);

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                for (int ei = 0; ei < 12; ei++)
                    ec(l,d,ijk)[ei] = hexa(vc(l,d,ijk)).edgeCenter(ei);
            }
        }
    }

    ec.correctAggData();

    return tEc;
}

// Instantiate

template class fvMeshMetrics<colocated>;
template class fvMeshMetrics<staggered>;

}

}

}
