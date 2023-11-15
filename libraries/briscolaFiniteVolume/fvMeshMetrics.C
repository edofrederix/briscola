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

    forAllLevels(fc, l, d, i, j, k)
    {
        const partLevelPoints& points = fvMsh_[l].points();

        // For each cell the face centers are calculated from the average of the
        // four face vertices (Wesseling, p. 483).

        const vector shift = MeshType::shift[d];
        const vector ijk(vector(i,j,k)+shift);

        fc(l,d,i,j,k).left() =
            0.25
          * (
                points.interp(ijk)
              + points.interp(ijk+vector(unitY))
              + points.interp(ijk+vector(unitZ))
              + points.interp(ijk+vector(unitYZ))
            );

        fc(l,d,i,j,k).right() =
            0.25
          * (
                points.interp(ijk+vector(unitX))
              + points.interp(ijk+vector(unitXY))
              + points.interp(ijk+vector(unitXZ))
              + points.interp(ijk+vector(unitXYZ))
            );

        fc(l,d,i,j,k).bottom() =
            0.25
          * (
                points.interp(ijk)
              + points.interp(ijk+vector(unitX))
              + points.interp(ijk+vector(unitZ))
              + points.interp(ijk+vector(unitXZ))
            );

        fc(l,d,i,j,k).top() =
            0.25
          * (
                points.interp(ijk+vector(unitY))
              + points.interp(ijk+vector(unitXY))
              + points.interp(ijk+vector(unitYZ))
              + points.interp(ijk+vector(unitXYZ))
            );

        fc(l,d,i,j,k).aft() =
            0.25
          * (
                points.interp(ijk)
              + points.interp(ijk+vector(unitX))
              + points.interp(ijk+vector(unitY))
              + points.interp(ijk+vector(unitXY))
            );

        fc(l,d,i,j,k).fore() =
            0.25
          * (
                points.interp(ijk+vector(unitZ))
              + points.interp(ijk+vector(unitXZ))
              + points.interp(ijk+vector(unitYZ))
              + points.interp(ijk+vector(unitXYZ))
            );
    }

    fc.correctParallelBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateEdgeCenters()
{
    meshField<edgeVector,MeshType>& ec = edgeCenters_;

    ec = Zero;

    forAllLevels(ec, l, d, i, j, k)
    {
        const partLevelPoints& points = fvMsh_[l].points();

        const vector shift = MeshType::shift[d];
        const vector ijk(vector(i,j,k)+shift);

        // Edges in x

        ec(l,d,i,j,k).ba() =
            0.5
          * (
                points.interp(ijk)
              + points.interp(ijk+vector(unitX))
            );

        ec(l,d,i,j,k).ta() =
            0.5
          * (
                points.interp(ijk+vector(unitY))
              + points.interp(ijk+vector(unitXY))
            );

        ec(l,d,i,j,k).bf() =
            0.5
          * (
                points.interp(ijk+vector(unitZ))
              + points.interp(ijk+vector(unitXZ))
            );

        ec(l,d,i,j,k).tf() =
            0.5
          * (
                points.interp(ijk+vector(unitYZ))
              + points.interp(ijk+vector(unitXYZ))
            );

        // Edges in y

        ec(l,d,i,j,k).la() =
            0.5
          * (
                points.interp(ijk)
              + points.interp(ijk+vector(unitY))
            );

        ec(l,d,i,j,k).ra() =
            0.5
          * (
                points.interp(ijk+vector(unitX))
              + points.interp(ijk+vector(unitXY))
            );

        ec(l,d,i,j,k).lf() =
            0.5
          * (
                points.interp(ijk+vector(unitZ))
              + points.interp(ijk+vector(unitYZ))
            );

        ec(l,d,i,j,k).rf() =
            0.5
          * (
                points.interp(ijk+vector(unitXZ))
              + points.interp(ijk+vector(unitXYZ))
            );

        // Edges in z

        ec(l,d,i,j,k).lb() =
            0.5
          * (
                points.interp(ijk)
              + points.interp(ijk+vector(unitZ))
            );

        ec(l,d,i,j,k).rb() =
            0.5
          * (
                points.interp(ijk+vector(unitX))
              + points.interp(ijk+vector(unitXZ))
            );

        ec(l,d,i,j,k).lt() =
            0.5
          * (
                points.interp(ijk+vector(unitY))
              + points.interp(ijk+vector(unitYZ))
            );

        ec(l,d,i,j,k).rt() =
            0.5
          * (
                points.interp(ijk+vector(unitXY))
              + points.interp(ijk+vector(unitXYZ))
            );
    }

    ec.correctParallelBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateVertexCenters()
{
    meshField<vertexVector,MeshType>& vc = vertexCenters_;

    vc = Zero;

    forAllLevels(vc, l, d, i, j, k)
    {
        const partLevelPoints& points = fvMsh_[l].points();

        const vector shift = MeshType::shift[d];
        const vector ijk(vector(i,j,k)+shift);

        vc(l,d,i,j,k).lba() = points.interp(ijk);
        vc(l,d,i,j,k).rba() = points.interp(ijk+vector(unitX));
        vc(l,d,i,j,k).lta() = points.interp(ijk+vector(unitY));
        vc(l,d,i,j,k).rta() = points.interp(ijk+vector(unitXY));
        vc(l,d,i,j,k).lbf() = points.interp(ijk+vector(unitZ));
        vc(l,d,i,j,k).rbf() = points.interp(ijk+vector(unitXZ));
        vc(l,d,i,j,k).ltf() = points.interp(ijk+vector(unitYZ));
        vc(l,d,i,j,k).rtf() = points.interp(ijk+vector(unitXYZ));
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

    forAllLevels(fn, l, d, i, j, k)
    {
        const partLevelPoints& points = fvMsh_[l].points();

        // For each cell the lower face normal in three directions is calculated
        // by taking half the cross product of the two vectors connecting the
        // diagonal vertex pairs (Wessling, p. 483). The normal's magnitude
        // equals the face area.

        const vector shift = MeshType::shift[d];

        const vector ijk(vector(i,j,k)+shift);

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

        fn(l,d,i,j,k) =
            faceVector
            (
                normalised(left),
                normalised(right),
                normalised(bottom),
                normalised(top),
                normalised(aft),
                normalised(fore)
            );

        fa(l,d,i,j,k) =
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

    forAllLevels(cc, l, d, i, j, k)
    {
        const partLevelPoints& points = fvMsh_[l].points();

        const vector shift = MeshType::shift[d];

        const vector ijk(vector(i,j,k)+shift);

        cc(l,d,i,j,k) =
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
            labelVector bo;

            for (bo.x() = -1; bo.x() <= 1; bo.x()++)
            for (bo.y() = -1; bo.y() <= 1; bo.y()++)
            for (bo.z() = -1; bo.z() <= 1; bo.z()++)
            if (cmptSum(cmptMag(bo)) > 0)
            {
                const labelVector S(cc.fvMsh().template S<MeshType>(l,d,bo));
                const labelVector E(cc.fvMsh().template E<MeshType>(l,d,bo));

                const label fi = faceNumber(bo);
                const label ei = edgeNumber(bo);
                const label vi = vertexNumber(bo);

                const label bod = cmptSum(cmptMag(bo));

                labelVector ijk;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    cc(l,d,ijk+bo) =
                        bod == 1 ? (2.0*fc(l,d,ijk)[fi] - cc(l,d,ijk))
                      : bod == 2 ? (2.0*ec(l,d,ijk)[ei] - cc(l,d,ijk))
                      :            (2.0*vc(l,d,ijk)[vi] - cc(l,d,ijk));
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
    const meshField<vertexVector,MeshType>& vc = vertexCenters_;

    // Small number to avoid division by zero in unset ghost cells

    cv = 1e-16;

    forAllLevels(cv, l, d, i, j, k)
    {
        if (fvMsh_.msh()[0].rectilinear() == unitXYZ)
        {
            // Cell volume is given by Wesseling's efficient formula (Wesseling,
            // p. 484)

            const vector Sx =
                fan(l,d,i,j,k).right() - fan(l,d,i,j,k).left();
            const vector Sy =
                fan(l,d,i,j,k).top()   - fan(l,d,i,j,k).bottom();
            const vector Sz =
                fan(l,d,i,j,k).fore()  - fan(l,d,i,j,k).aft();

            const vector Dx =
                fc(l,d,i,j,k).right() - fc(l,d,i,j,k).left();
            const vector Dy =
                fc(l,d,i,j,k).top()   - fc(l,d,i,j,k).bottom();
            const vector Dz =
                fc(l,d,i,j,k).fore()  - fc(l,d,i,j,k).aft();

            cv(l,d,i,j,k) = ((Dx & Sx) + (Dy & Sy) + (Dz & Sz))/6.0;
        }
        else
        {
            // Cell volume computed from tet decomposition, needed for
            // consistency with geometric VoF calculations.

            const vertexVector& v = vc(l,d,i,j,k);

            cv(l,d,i,j,k) = 0.0;

            for (int t = 0; t < numberOfTets; t++)
                cv(l,d,i,j,k) +=
                    tetVolume
                    (
                        v[tetDecomp[t][0]],
                        v[tetDecomp[t][1]],
                        v[tetDecomp[t][2]],
                        v[tetDecomp[t][3]]
                    );
        }
    }

    // Extrapolate cell volumes to ghost cells

    forAll(cv, l)
    {
        const labelVector N(fvMsh_[l].N());

        forAll(cv[l], d)
        {
            labelVector bo;

            for (bo.x() = -1; bo.x() <= 1; bo.x()++)
            for (bo.y() = -1; bo.y() <= 1; bo.y()++)
            for (bo.z() = -1; bo.z() <= 1; bo.z()++)
            if (cmptSum(cmptMag(bo)) > 0)
            {
                const labelVector S(cv.fvMsh().template S<MeshType>(l,d,bo));
                const labelVector E(cv.fvMsh().template E<MeshType>(l,d,bo));

                labelVector ijk;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    cv(l,d,ijk+bo) = cv(l,d,ijk);
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

    forAllLevels(cc, l, d, i, j, k)
    {
        fd(l,d,i,j,k) =
            faceScalar
            (
                1.0/Foam::mag(cc(l,d,i,j,k)-cc(l,d,i-1,j,k)),
                1.0/Foam::mag(cc(l,d,i,j,k)-cc(l,d,i+1,j,k)),
                1.0/Foam::mag(cc(l,d,i,j,k)-cc(l,d,i,j-1,k)),
                1.0/Foam::mag(cc(l,d,i,j,k)-cc(l,d,i,j+1,k)),
                1.0/Foam::mag(cc(l,d,i,j,k)-cc(l,d,i,j,k-1)),
                1.0/Foam::mag(cc(l,d,i,j,k)-cc(l,d,i,j,k+1))
            );
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
