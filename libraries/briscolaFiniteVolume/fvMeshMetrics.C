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
void fvMeshMetrics<MeshType>::calculateFaceCenters()
{
    meshField<faceVector,MeshType>& fc = faceCenters_;

    fc = Zero;

    forAll(fvMsh_, l)
    {
        const partLevelPoints& points = fvMsh_[l].points();

        // For each cell the face centers are calculated from the average of the
        // four face vertices (Wesseling, p. 483).

        forAll(fc[l], d)
        {
            meshDirection<faceVector,MeshType>& fcld = fc[l][d];

            const vector shift = MeshType::shift[d];

            forAllCells(fcld, i, j, k)
            {
                vector ijk(vector(i,j,k)+shift);

                fcld(i,j,k).left() =
                    0.25
                  * (
                        points.interp(ijk)
                      + points.interp(ijk+vector(unitY))
                      + points.interp(ijk+vector(unitZ))
                      + points.interp(ijk+vector(unitYZ))
                    );

                fcld(i,j,k).right() =
                    0.25
                  * (
                        points.interp(ijk+vector(unitX))
                      + points.interp(ijk+vector(unitXY))
                      + points.interp(ijk+vector(unitXZ))
                      + points.interp(ijk+vector(unitXYZ))
                    );

                fcld(i,j,k).bottom() =
                    0.25
                  * (
                        points.interp(ijk)
                      + points.interp(ijk+vector(unitX))
                      + points.interp(ijk+vector(unitZ))
                      + points.interp(ijk+vector(unitXZ))
                    );

                fcld(i,j,k).top() =
                    0.25
                  * (
                        points.interp(ijk+vector(unitY))
                      + points.interp(ijk+vector(unitXY))
                      + points.interp(ijk+vector(unitYZ))
                      + points.interp(ijk+vector(unitXYZ))
                    );

                fcld(i,j,k).aft() =
                    0.25
                  * (
                        points.interp(ijk)
                      + points.interp(ijk+vector(unitX))
                      + points.interp(ijk+vector(unitY))
                      + points.interp(ijk+vector(unitXY))
                    );

                fcld(i,j,k).fore() =
                    0.25
                  * (
                        points.interp(ijk+vector(unitZ))
                      + points.interp(ijk+vector(unitXZ))
                      + points.interp(ijk+vector(unitYZ))
                      + points.interp(ijk+vector(unitXYZ))
                    );
            }
        }
    }

    fc.correctParallelBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateEdgeCenters()
{
    meshField<edgeVector,MeshType>& ec = edgeCenters_;

    ec = Zero;

    forAll(fvMsh_, l)
    {
        const partLevelPoints& points = fvMsh_[l].points();

        forAll(ec[l], d)
        {
            meshDirection<edgeVector,MeshType>& ecld = ec[l][d];

            const vector shift = MeshType::shift[d];

            forAllCells(ecld, i, j, k)
            {
                vector ijk(vector(i,j,k)+shift);

                // Edges in x

                ecld(i,j,k).ba() =
                    0.5
                  * (
                        points.interp(ijk)
                      + points.interp(ijk+vector(unitX))
                    );

                ecld(i,j,k).ta() =
                    0.5
                  * (
                        points.interp(ijk+vector(unitY))
                      + points.interp(ijk+vector(unitXY))
                    );

                ecld(i,j,k).bf() =
                    0.5
                  * (
                        points.interp(ijk+vector(unitZ))
                      + points.interp(ijk+vector(unitXZ))
                    );

                ecld(i,j,k).tf() =
                    0.5
                  * (
                        points.interp(ijk+vector(unitYZ))
                      + points.interp(ijk+vector(unitXYZ))
                    );

                // Edges in y

                ecld(i,j,k).la() =
                    0.5
                  * (
                        points.interp(ijk)
                      + points.interp(ijk+vector(unitY))
                    );

                ecld(i,j,k).ra() =
                    0.5
                  * (
                        points.interp(ijk+vector(unitX))
                      + points.interp(ijk+vector(unitXY))
                    );

                ecld(i,j,k).lf() =
                    0.5
                  * (
                        points.interp(ijk+vector(unitZ))
                      + points.interp(ijk+vector(unitYZ))
                    );

                ecld(i,j,k).rf() =
                    0.5
                  * (
                        points.interp(ijk+vector(unitXZ))
                      + points.interp(ijk+vector(unitXYZ))
                    );

                // Edges in z

                ecld(i,j,k).lb() =
                    0.5
                  * (
                        points.interp(ijk)
                      + points.interp(ijk+vector(unitZ))
                    );

                ecld(i,j,k).rb() =
                    0.5
                  * (
                        points.interp(ijk+vector(unitX))
                      + points.interp(ijk+vector(unitXZ))
                    );

                ecld(i,j,k).lt() =
                    0.5
                  * (
                        points.interp(ijk+vector(unitY))
                      + points.interp(ijk+vector(unitYZ))
                    );

                ecld(i,j,k).rt() =
                    0.5
                  * (
                        points.interp(ijk+vector(unitXY))
                      + points.interp(ijk+vector(unitXYZ))
                    );
            }
        }
    }

    ec.correctParallelBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateVertexCenters()
{
    meshField<vertexVector,MeshType>& vc = vertexCenters_;

    vc = Zero;

    forAll(fvMsh_, l)
    {
        const partLevelPoints& points = fvMsh_[l].points();

        forAll(vc[l], d)
        {
            meshDirection<vertexVector,MeshType>& vcld = vc[l][d];

            const vector shift = MeshType::shift[d];

            forAllCells(vcld, i, j, k)
            {
                vector ijk(vector(i,j,k)+shift);

                vcld(i,j,k).lba() = points.interp(ijk);
                vcld(i,j,k).rba() = points.interp(ijk+vector(unitX));
                vcld(i,j,k).lta() = points.interp(ijk+vector(unitY));
                vcld(i,j,k).rta() = points.interp(ijk+vector(unitXY));
                vcld(i,j,k).lbf() = points.interp(ijk+vector(unitZ));
                vcld(i,j,k).rbf() = points.interp(ijk+vector(unitXZ));
                vcld(i,j,k).ltf() = points.interp(ijk+vector(unitYZ));
                vcld(i,j,k).rtf() = points.interp(ijk+vector(unitXYZ));
            }
        }
    }

    vc.correctParallelBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceAreasAndNormals()
{
    meshField<faceVector,MeshType>& fn = faceNormals_;
    meshField<faceScalar,MeshType>& fa = faceAreas_;
    meshField<faceVector,MeshType>& fan = faceAreaNormals_;

    fn = Zero;
    fa = Zero;
    fan = Zero;

    forAll(fn, l)
    {
        const partLevelPoints& points = fvMsh_[l].points();

        // For each cell the lower face normal in three directions is calculated
        // by taking half the cross product of the two vectors connecting the
        // diagonal vertex pairs (Wessling, p. 483). The normal's magnitude
        // equals the face area.

        forAll(fn[l], d)
        {
            meshDirection<faceVector,MeshType>& fnld = fn[l][d];
            meshDirection<faceScalar,MeshType>& fald = fa[l][d];

            const vector shift = MeshType::shift[d];

            forAllCells(fnld, i, j, k)
            {
                vector ijk(vector(i,j,k)+shift);

                const vector left =
                  - 0.5
                  * (
                        (
                            points.interp(ijk+vector(unitYZ))
                          - points.interp(ijk)
                        )
                      ^ (
                            points.interp(ijk+vector(unitZ))
                          - points.interp(ijk+vector(unitY))
                        )
                    );

                const vector right =
                    0.5
                  * (
                        (
                            points.interp(ijk+vector(unitXYZ))
                          - points.interp(ijk+vector(unitX))
                        )
                      ^ (
                            points.interp(ijk+vector(unitXZ))
                          - points.interp(ijk+vector(unitXY))
                        )
                    );

                const vector bottom =
                  - 0.5
                  * (
                        (
                            points.interp(ijk+vector(unitZ))
                          - points.interp(ijk+vector(unitX))
                        )
                      ^ (
                            points.interp(ijk+vector(unitXZ))
                          - points.interp(ijk)
                        )
                    );

                const vector top =
                    0.5
                  * (
                        (
                            points.interp(ijk+vector(unitYZ))
                          - points.interp(ijk+vector(unitXY))
                        )
                      ^ (
                            points.interp(ijk+vector(unitXYZ))
                          - points.interp(ijk+vector(unitY))
                        )
                    );

                const vector aft =
                  - 0.5
                  * (
                        (
                            points.interp(ijk+vector(unitXY))
                          - points.interp(ijk)
                        )
                      ^ (
                            points.interp(ijk+vector(unitY))
                          - points.interp(ijk+vector(unitX))
                        )
                    );

                const vector fore =
                    0.5
                  * (
                        (
                            points.interp(ijk+vector(unitXYZ))
                          - points.interp(ijk+vector(unitZ))
                        )
                      ^ (
                            points.interp(ijk+vector(unitYZ))
                          - points.interp(ijk+vector(unitXZ))
                        )
                    );

                fnld(i,j,k) =
                    faceVector
                    (
                        normalised(left),
                        normalised(right),
                        normalised(bottom),
                        normalised(top),
                        normalised(aft),
                        normalised(fore)
                    );

                fald(i,j,k) =
                    faceScalar
                    (
                        mag(left),
                        mag(right),
                        mag(bottom),
                        mag(top),
                        mag(aft),
                        mag(fore)
                    );
            }
        }
    }

    fa.correctCommBoundaryConditions();
    fn.correctCommBoundaryConditions();

    fan = fa*fn;

    fan.correctCommBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateCellCenters()
{
    meshField<vector,MeshType>& cc = cellCenters_;

    cc = Zero;

    const meshField<faceVector,MeshType>& fc = faceCenters_;
    const meshField<edgeVector,MeshType>& ec = edgeCenters_;
    const meshField<vertexVector,MeshType>& vc = vertexCenters_;

    // First, set internal cell centers from the point coordinates

    forAll(cc, l)
    {
        const partLevelPoints& points = fvMsh_[l].points();

        forAll(cc[l], d)
        {
            meshDirection<vector,MeshType>& ccld = cc[l][d];

            const vector shift = MeshType::shift[d];

            forAllCells(ccld, i, j, k)
            {
                vector ijk(vector(i,j,k)+shift);

                ccld(i,j,k) =
                    0.125
                  * (
                        points.interp(ijk)
                      + points.interp(ijk+vector(unitX))
                      + points.interp(ijk+vector(unitY))
                      + points.interp(ijk+vector(unitXY))
                      + points.interp(ijk+vector(unitZ))
                      + points.interp(ijk+vector(unitXZ))
                      + points.interp(ijk+vector(unitYZ))
                      + points.interp(ijk+vector(unitXYZ))
                    );
            }
        }
    }

    // Correct ghost cells. We have the following situations that need to be
    // handled:
    //
    //  1) All normal (non-communicating) boundaries have their cell centers
    //     computed by a cell-center-to-point vector projection, which is
    //     performed first for all ghost cells. In principle some ghost cell
    //     centers (especially for colocated grids) could be computed from the
    //     points directly, however, these points are also projections so this
    //     will give the same result.
    //
    //  2) Parallel boundaries are set by correctParallelBoundaryConditions().
    //
    //  3) Periodic boundaries are set using the same projection as for normal
    //     boundaries. This is only accurate for symmetry across the periodic
    //     boundary, and could be improved by communicating cell-center-to-point
    //     vectors.

    // First set all ghost cells using a cell-center-to-point projection

    forAll(cc, l)
    {
        const labelVector N(fvMsh_[l].N());

        forAll(cc[l], d)
        {
            meshDirection<vector,MeshType>& ccld = cc[l][d];

            labelVector bo;

            for (bo.x() = -1; bo.x() <= 1; bo.x()++)
            for (bo.y() = -1; bo.y() <= 1; bo.y()++)
            for (bo.z() = -1; bo.z() <= 1; bo.z()++)
            if (cmptSum(cmptMag(bo)) > 0)
            {
                const labelVector S(ccld.boundaryStart(bo));
                const labelVector E(ccld.boundaryEnd(bo));

                const label fi = faceNumber(bo);
                const label ei = edgeNumber(bo);
                const label vi = vertexNumber(bo);

                const label bod = cmptSum(cmptMag(bo));

                labelVector ijk;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    ccld(ijk+bo) =
                        bod == 1 ? (2.0*fc[l][d](ijk)[fi] - ccld(ijk))
                      : bod == 2 ? (2.0*ec[l][d](ijk)[ei] - ccld(ijk))
                      :            (2.0*vc[l][d](ijk)[vi] - ccld(ijk));
                }
            }
        }
    }

    // Overwrite with a parallel update

    cc.correctParallelBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateCellVolumes()
{
    meshField<scalar,MeshType>& cv = cellVolumes_;

    const meshField<faceVector,MeshType>& fan = faceAreaNormals_;
    const meshField<faceVector,MeshType>& fc = faceCenters_;

    // Small number to avoid division by zero in unset ghost cells

    cv = 1e-16;

    forAll(cv, l)
    {
        // Cell volume is given by Wesseling's efficient formula (Wesseling, p.
        // 484)

        forAll(cv[l], d)
        {
            meshDirection<scalar,MeshType>& cvld = cv[l][d];

            const meshDirection<faceVector,MeshType>& fanld = fan[l][d];
            const meshDirection<faceVector,MeshType>& fcld = fc[l][d];

            forAllCells(cvld, i, j, k)
            {
                const vector Sx = fanld(i,j,k).right() - fanld(i,j,k).left();
                const vector Sy = fanld(i,j,k).top()   - fanld(i,j,k).bottom();
                const vector Sz = fanld(i,j,k).fore()  - fanld(i,j,k).aft();

                const vector Dx = fcld(i,j,k).right() - fcld(i,j,k).left();
                const vector Dy = fcld(i,j,k).top()   - fcld(i,j,k).bottom();
                const vector Dz = fcld(i,j,k).fore()  - fcld(i,j,k).aft();

                cvld(i,j,k) = ((Dx & Sx) + (Dy & Sy) + (Dz & Sz))/6.0;
            }
        }
    }

    // Extrapolate cell volumes to ghost cells

    forAll(cv, l)
    {
        const labelVector N(fvMsh_[l].N());

        forAll(cv[l], d)
        {
            meshDirection<scalar,MeshType>& cvld = cv[l][d];

            labelVector bo;

            for (bo.x() = -1; bo.x() <= 1; bo.x()++)
            for (bo.y() = -1; bo.y() <= 1; bo.y()++)
            for (bo.z() = -1; bo.z() <= 1; bo.z()++)
            if (cmptSum(cmptMag(bo)) > 0)
            {
                const labelVector S(cvld.boundaryStart(bo));
                const labelVector E(cvld.boundaryEnd(bo));

                labelVector ijk;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    cvld(ijk+bo) = cvld(ijk);
                }
            }
        }
    }

    cv.correctCommBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceDeltas()
{
    meshField<faceScalar,MeshType>& fd = faceDeltas_;

    const meshField<vector,MeshType>& cc = cellCenters_;

    fd = Zero;

    forAll(cc, l)
    {
        forAll(cc[l], d)
        {
            meshDirection<faceScalar,MeshType>& fdld = fd[l][d];
            const meshDirection<vector,MeshType>& ccld = cc[l][d];

            forAllCells(fdld, i, j, k)
            {
                fdld(i,j,k) =
                    faceScalar
                    (
                        1.0/Foam::mag(ccld(i,j,k)-ccld(i-1,j,k)),
                        1.0/Foam::mag(ccld(i,j,k)-ccld(i+1,j,k)),
                        1.0/Foam::mag(ccld(i,j,k)-ccld(i,j-1,k)),
                        1.0/Foam::mag(ccld(i,j,k)-ccld(i,j+1,k)),
                        1.0/Foam::mag(ccld(i,j,k)-ccld(i,j,k-1)),
                        1.0/Foam::mag(ccld(i,j,k)-ccld(i,j,k+1))
                    );
            }
        }
    }

    fd.correctCommBoundaryConditions();
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
        word(MeshType::typeName) + "FaceCenters",
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
    )
{
    calculateFaceCenters();
    calculateEdgeCenters();
    calculateVertexCenters();
    calculateFaceAreasAndNormals();
    calculateCellCenters();
    calculateCellVolumes();
    calculateFaceDeltas();
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
