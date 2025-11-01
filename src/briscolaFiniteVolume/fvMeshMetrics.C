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
void fvMeshMetrics<MeshType>::calculateVertexCenters()
{
    meshField<vertexVector,MeshType>& vc = vertexCenters_;

    vc = Zero;

    forAll(fvMsh_, l)
    {
        const partPoints& points = fvMsh_[l].points();

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);
            const vector shift = MeshType::shift[d];

            // Create vertices from a partPoints object, which consistently
            // generates ghost points

            partPoints p(fvMsh_.msh());

            p.clear();
            p.setSizeFromCells(N);

            forAllBlock(p, i, j, k)
                p(i,j,k) = points.interp(vector(i,j,k) + shift);

            p.calcGhostPoints();
            p.clean();

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
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceCenters()
{
    const meshField<faceVector,MeshType> fcAoS(aos().faceCenters());

    faceField<vector,MeshType>& fc = faceCenters_;

    fc = Zero;

    forAllFaces(fc, fd, l, d, i, j, k)
        fc[fd](l,d,i,j,k) = fcAoS(l,d,i,j,k)[fd*2];
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceAreasAndNormals()
{
    const meshField<faceVector,MeshType> fnAoS(aos().faceNormals());
    const meshField<faceScalar,MeshType> faAoS(aos().faceAreas());
    const meshField<faceVector,MeshType> fanAoS(aos().faceAreaNormals());

    faceField<vector,MeshType>& fn = faceNormals_;
    faceField<scalar,MeshType>& fa = faceAreas_;
    faceField<vector,MeshType>& fan = faceAreaNormals_;

    fn = Zero;
    fa = Zero;
    fan = Zero;

    forAllFaces(fn, fd, l, d, i, j, k)
    {
        fn[fd](l,d,i,j,k) = fnAoS(l,d,i,j,k)[fd*2];
        fa[fd](l,d,i,j,k) = faAoS(l,d,i,j,k)[fd*2];
        fan[fd](l,d,i,j,k) = fanAoS(l,d,i,j,k)[fd*2];
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateCellCenters()
{
    const meshField<vertexVector,MeshType>& vc = vertexCenters_;

    meshField<vector,MeshType>& cc = cellCenters_;

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
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateCellVolumes()
{
    const meshField<vertexVector,MeshType>& vc = vertexCenters_;

    meshField<scalar,MeshType>& cv = cellVolumes_;

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

    // Compute at store the inverse cell volumes
    inverseCellVolumes_ = 1.0/cv;
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceDeltas()
{
    const meshField<faceScalar,MeshType> deltaAoS(aos().faceDeltas());

    faceField<scalar,MeshType>& delta = faceDeltas_;

    delta = Zero;

    forAllFaces(delta, fd, l, d, i, j, k)
    {
        delta[fd](l,d,i,j,k) = deltaAoS(l,d,i,j,k)[fd*2];
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceWeights()
{
    const meshField<faceScalar,MeshType> fwcAoS(aos().faceWeightsCenter());
    const meshField<faceScalar,MeshType> fwnAoS(aos().faceWeightsNeighbor());

    faceField<scalar,MeshType>& fwc = faceWeightsCenter_;
    faceField<scalar,MeshType>& fwn = faceWeightsNeighbor_;

    // Initialize weights as 50/50

    fwc = 0.5;
    fwn = 0.5;

    forAllFaces(fwc, fd, l, d, i, j, k)
    {
        fwc[fd](l,d,i,j,k) = fwcAoS(l,d,i,j,k)[fd*2];
        fwn[fd](l,d,i,j,k) = fwnAoS(l,d,i,j,k)[fd*2];
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setGlobalCellNumbers()
{
    globalCellNumbers_ = -1;
    globalCellCounts_.resize(globalCellNumbers_.size());

    // Cells are numbered lexicographically on each processor and with
    // increasing processor rank

    forAll(globalCellNumbers_, l)
    {
        globalCellCounts_[l].resize(MeshType::numberOfDirections);

        forAll(globalCellNumbers_[l], d)
        {
            labelList sizes(Pstream::nProcs());

            sizes[Pstream::myProcNo()] =
                cmptProduct(fvMsh_.N<MeshType>(l,d));

            Pstream::gatherList(sizes);
            Pstream::scatterList(sizes);

            label start = 0;

            for (int i = 0; i < Pstream::myProcNo(); i++)
                start += sizes[i];

            int c = 0;
            forAllCells(globalCellNumbers_[l][d], i, j, k)
            {
                globalCellNumbers_[l][d](i,j,k) = start + c++;
            }

            globalCellCounts_[l][d] = sum(sizes);
        }
    }

    globalCellNumbers_.correctCommsBoundaryConditions();
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
    mask_ = Zero;

    if (immersedBoundaries_.size())
    {
        forAll(immersedBoundaries_, i)
            mask_ += immersedBoundaries_[i].mask();

        mask_ = min(mask_, 1.0);
    }
}

template<class MeshType>
fvMeshMetrics<MeshType>::fvMeshMetrics(const fvMesh& fvMsh)
:
    fvMsh_(fvMsh),
    faceCenters_
    (
        word(MeshType::typeName) + "FaceCenters",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    vertexCenters_
    (
        word(MeshType::typeName) + "VertexCenters",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    faceNormals_
    (
        word(MeshType::typeName) + "FaceNormals",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    faceAreas_
    (
        word(MeshType::typeName) + "FaceAreas",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    faceAreaNormals_
    (
        word(MeshType::typeName) + "FaceAreaNormals",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    cellCenters_
    (
        word(MeshType::typeName) + "CellCenters",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    cellVolumes_
    (
        word(MeshType::typeName) + "CellVolumes",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    inverseCellVolumes_
    (
        word(MeshType::typeName) + "InverseCellVolumes",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    faceDeltas_
    (
        word(MeshType::typeName) + "FaceDeltas",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    faceWeightsCenter_
    (
        word(MeshType::typeName) + "FaceWeightsCenter",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    faceWeightsNeighbor_
    (
        word(MeshType::typeName) + "FaceWeightsNeighbor",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    globalCellNumbers_
    (
        word(MeshType::typeName) + "GlobalCellNumbers",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    immersedBoundaries_(),
    mask_
    (
        word(MeshType::typeName) + "Mask",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    aosPtr_(new AoS(*this))
{
    calculateVertexCenters();
    calculateFaceCenters();
    calculateFaceAreasAndNormals();
    calculateCellCenters();
    calculateCellVolumes();
    calculateFaceDeltas();
    calculateFaceWeights();
    setGlobalCellNumbers();
    setImmersedBoundaries();
    setImmersedBoundaryMask();
}

// On-demand Array of Structures (AoS) metrics class

template<class MeshType>
tmp<meshField<faceVector,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceCenters() const
{
    const meshField<vertexVector,MeshType>& vc = metrics_.vertexCenters_;

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

    return tFc;
}

template<class MeshType>
tmp<meshField<edgeVector,MeshType>>
fvMeshMetrics<MeshType>::AoS::edgeCenters() const
{
    const meshField<vertexVector,MeshType>& vc = metrics_.vertexCenters_;

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

    return tEc;
}

template<class MeshType>
tmp<meshField<faceVector,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceNormals() const
{
    const meshField<vertexVector,MeshType>& vc = metrics_.vertexCenters_;

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

    return tFn;
}

template<class MeshType>
tmp<meshField<faceScalar,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceAreas() const
{
    const meshField<vertexVector,MeshType>& vc = metrics_.vertexCenters_;

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

    return tFa;
}

template<class MeshType>
tmp<meshField<faceVector,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceAreaNormals() const
{
    const meshField<vertexVector,MeshType>& vc = metrics_.vertexCenters_;

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

    return tFan;
}

template<class MeshType>
tmp<meshField<faceScalar,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceDeltas() const
{
    const meshField<vector,MeshType>& cc = metrics_.cellCenters_;

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

    return tDelta;
}

template<class MeshType>
tmp<meshField<faceScalar,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceWeightsCenter() const
{
    const meshField<vector,MeshType>& cc = metrics_.cellCenters_;

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

    return tFwc;
}

template<class MeshType>
tmp<meshField<faceScalar,MeshType>>
fvMeshMetrics<MeshType>::AoS::faceWeightsNeighbor() const
{
    const meshField<vector,MeshType>& cc = metrics_.cellCenters_;

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

    return tFwn;
}

// Instantiate

template class fvMeshMetrics<colocated>;
template class fvMeshMetrics<staggered>;

}

}

}
