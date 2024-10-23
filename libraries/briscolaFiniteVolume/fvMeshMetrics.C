#include "fvMeshMetrics.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

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
            // generated ghost points

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
            }
        }
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceCenters()
{
    const meshField<vertexVector,MeshType>& vc = vertexCenters_;

    meshField<faceVector,MeshType>& fc = faceCenters_;

    fc = Zero;

    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            // For each cell the face centers are calculated from the average of
            // the four face vertices (Wesseling, p. 483).

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                for (int fi = 0; fi < 6; fi++)
                    for (int vi = 0; vi < 4; vi++)
                        fc(l,d,ijk)[fi] +=
                            0.25*vc(l,d,ijk)[vertexNumsInFace[fi][vi]];
            }
        }
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateEdgeCenters()
{
    const meshField<vertexVector,MeshType>& vc = vertexCenters_;

    meshField<edgeVector,MeshType>& ec = edgeCenters_;

    ec = Zero;

    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            // For each cell the edge centers are calculated from the average of
            // the two edge vertices

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                for (int ei = 0; ei < 12; ei++)
                    for (int vi = 0; vi < 2; vi++)
                        ec(l,d,ijk)[ei] +=
                            0.5*vc(l,d,ijk)[vertexNumsInEdge[ei][vi]];
            }
        }
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceAreasAndNormals()
{
    const meshField<vertexVector,MeshType>& vc = vertexCenters_;

    meshField<faceVector,MeshType>& fn = faceNormals_;
    meshField<faceScalar,MeshType>& fa = faceAreas_;
    meshField<faceVector,MeshType>& fan = faceAreaNormals_;

    fn = Zero;
    fa = Zero;
    fan = Zero;

    forAll(fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector N = fvMsh_.N<MeshType>(l,d);

            // For each cell the outward face normals are calculated by taking
            // half the cross product of the two vectors connecting the diagonal
            // vertex pairs (Wesseling, p. 483). The normal's magnitude equals
            // the face area.

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                for (int fi = 0; fi < 6; fi++)
                {
                    fan(l,d,ijk)[fi] =
                        0.5
                      * (
                            (
                                vc(l,d,ijk)[vertexNumsInFaceCC[fi][2]]
                              - vc(l,d,ijk)[vertexNumsInFaceCC[fi][0]]
                            )
                          ^ (
                                vc(l,d,ijk)[vertexNumsInFaceCC[fi][3]]
                              - vc(l,d,ijk)[vertexNumsInFaceCC[fi][1]]
                            )
                        );

                    fa(l,d,ijk)[fi] = Foam::mag(fan(l,d,ijk)[fi]);
                    fn(l,d,ijk)[fi] = Foam::normalised(fan(l,d,ijk)[fi]);
                }
            }
        }
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

            // For each cell the cell center is calculated from the average of
            // the eight vertices

            labelVector ijk;
            for (ijk.x() = -1; ijk.x() < N.x() + 1; ijk.x()++)
            for (ijk.y() = -1; ijk.y() < N.y() + 1; ijk.y()++)
            for (ijk.z() = -1; ijk.z() < N.z() + 1; ijk.z()++)
            {
                for (int vi = 0; vi < 8; vi++)
                    cc(l,d,ijk) += 0.125*vc(l,d,ijk)[vi];
            }
        }
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateCellVolumes()
{
    const meshField<vertexVector,MeshType>& vc = vertexCenters_;

    meshField<scalar,MeshType>& cv = cellVolumes_;

    // Small number to avoid division by zero in unset ghost cells

    cv = 1e-16;

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
                for (int ti = 0; ti < numberOfTets; ti++)
                    cv(l,d,ijk) +=
                        tetVolume
                        (
                            vc(l,d,ijk)[tetDecomp[ti][0]],
                            vc(l,d,ijk)[tetDecomp[ti][1]],
                            vc(l,d,ijk)[tetDecomp[ti][2]],
                            vc(l,d,ijk)[tetDecomp[ti][3]]
                        );
            }
        }
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceDeltas()
{
    meshField<faceScalar,MeshType>& delta = faceDeltas_;

    const meshField<vector,MeshType>& cc = cellCenters_;

    delta = Zero;

    forAllFaces(delta, l, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNei(i,j,k,fd));

        delta(l,d,ijk)[fd*2] =
            1.0/Foam::mag(cc(l,d,ijk)-cc(l,d,nei));

        delta(l,d,nei)[fd*2+1] = delta(l,d,ijk)[fd*2];
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceWeights()
{
    meshField<faceScalar,MeshType>& fwc = faceWeightsCenter_;
    meshField<faceScalar,MeshType>& fwn = faceWeightsNeighbor_;

    const meshField<faceScalar,MeshType>& delta = faceDeltas_;
    const meshField<vector,MeshType>& cc = cellCenters_;
    const meshField<faceVector,MeshType>& fc = faceCenters_;

    fwc = Zero;
    fwn = Zero;

    forAllFaces(fwc, l, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNei(i,j,k,fd));

        fwc(l,d,ijk)[fd*2] =
            Foam::mag(cc(l,d,nei) - fc(l,d,ijk)[fd*2])*delta(l,d,ijk)[fd*2];

        fwn(l,d,ijk)[fd*2] =
            Foam::mag(cc(l,d,ijk) - fc(l,d,ijk)[fd*2])*delta(l,d,ijk)[fd*2];

        fwc(l,d,nei)[fd*2+1] = fwn(l,d,ijk)[fd*2];
        fwn(l,d,nei)[fd*2+1] = fwc(l,d,ijk)[fd*2];
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::setGlobalCellNumbers()
{
    globalCellNumbers_ = -1;

    forAll(globalCellNumbers_, l)
    {
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
        }
    }

    globalCellNumbers_.correctCommsBoundaryConditions();
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
    edgeCenters_
    (
        word(MeshType::typeName) + "EdgeCenters",
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
    )
{
    calculateVertexCenters();
    calculateFaceCenters();
    calculateEdgeCenters();
    calculateFaceAreasAndNormals();
    calculateCellCenters();
    calculateCellVolumes();
    calculateFaceDeltas();
    calculateFaceWeights();
    setGlobalCellNumbers();
}

template<class MeshType>
fvMeshMetrics<MeshType>::~fvMeshMetrics()
{}

defineTemplateTypeNameAndDebug(fvMeshMetrics<colocated>, 0);
defineTemplateTypeNameAndDebug(fvMeshMetrics<staggered>, 0);

template class fvMeshMetrics<colocated>;
template class fvMeshMetrics<staggered>;

}

}

}
