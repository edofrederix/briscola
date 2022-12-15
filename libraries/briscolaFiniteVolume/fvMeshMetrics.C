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

    forAll(fvMsh_, leveli)
    {
        const partLevelPoints& points = fvMsh_[leveli].points();

        // For each cell the face centers are calculated from the average of the
        // four face vertices (Wesseling, p. 483).

        forAll(fc[leveli], d)
        {
            meshDirection<faceVector,MeshType>& fcld = fc[leveli][d];

            const vector shift = MeshType::shift[d];

            forAllBlock(fcld, i, j, k)
            {
                fcld(i,j,k).left() =
                    0.25
                  * (
                        points.interp(i,j,  k,  shift)
                      + points.interp(i,j+1,k,  shift)
                      + points.interp(i,j,  k+1,shift)
                      + points.interp(i,j+1,k+1,shift)
                    );

                fcld(i,j,k).right() =
                    0.25
                  * (
                        points.interp(i+1,j,  k,  shift)
                      + points.interp(i+1,j+1,k,  shift)
                      + points.interp(i+1,j,  k+1,shift)
                      + points.interp(i+1,j+1,k+1,shift)
                    );

                fcld(i,j,k).bottom() =
                    0.25
                  * (
                        points.interp(i,  j,k,  shift)
                      + points.interp(i+1,j,k,  shift)
                      + points.interp(i,  j,k+1,shift)
                      + points.interp(i+1,j,k+1,shift)
                    );

                fcld(i,j,k).top() =
                    0.25
                  * (
                        points.interp(i,  j+1,k,  shift)
                      + points.interp(i+1,j+1,k,  shift)
                      + points.interp(i,  j+1,k+1,shift)
                      + points.interp(i+1,j+1,k+1,shift)
                    );

                fcld(i,j,k).aft() =
                    0.25
                  * (
                        points.interp(i,  j,  k,shift)
                      + points.interp(i+1,j,  k,shift)
                      + points.interp(i,  j+1,k,shift)
                      + points.interp(i+1,j+1,k,shift)
                    );

                fcld(i,j,k).fore() =
                    0.25
                  * (
                        points.interp(i,  j,  k+1,shift)
                      + points.interp(i+1,j,  k+1,shift)
                      + points.interp(i,  j+1,k+1,shift)
                      + points.interp(i+1,j+1,k+1,shift)
                    );
            }
        }
    }

    fc.correctBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceAreasAndNormals()
{
    meshField<faceVector,MeshType>& fn = faceNormals_;
    meshField<faceScalar,MeshType>& fa = faceAreas_;
    meshField<faceVector,MeshType>& fan = faceAreaNormals_;

    forAll(fn, leveli)
    {
        const partLevelPoints& points = fvMsh_[leveli].points();

        // For each cell the lower face normal in three directions is calculated
        // by taking half the cross product of the two vectors connecting the
        // diagonal vertex pairs (Wessling, p. 483). The normal's magnitude
        // equals the face area.

        forAll(fn[leveli], d)
        {
            meshDirection<faceVector,MeshType>& fnld = fn[leveli][d];
            meshDirection<faceScalar,MeshType>& fald = fa[leveli][d];

            const vector shift = MeshType::shift[d];

            forAllBlock(fnld, i, j, k)
            {
                const vector left =
                    0.5
                  * (
                        (
                            points.interp(i,j+1,k+1,shift)
                          - points.interp(i,j,  k,  shift)
                        )
                      ^ (
                            points.interp(i,j,  k+1,shift)
                          - points.interp(i,j+1,k,  shift)
                        )
                    );

                const vector right =
                    0.5
                  * (
                        (
                            points.interp(i+1,j+1,k+1,shift)
                          - points.interp(i+1,j,  k,  shift)
                        )
                      ^ (
                            points.interp(i+1,j,  k+1,shift)
                          - points.interp(i+1,j+1,k,  shift)
                        )
                    );

                const vector bottom =
                    0.5
                  * (
                        (
                            points.interp(i,  j,k+1,shift)
                          - points.interp(i+1,j,k,  shift)
                        )
                      ^ (
                            points.interp(i+1,j,k+1,shift)
                          - points.interp(i,  j,k,  shift)
                        )
                    );

                const vector top =
                    0.5
                  * (
                        (
                            points.interp(i,  j+1,k+1,shift)
                          - points.interp(i+1,j+1,k,  shift)
                        )
                      ^ (
                            points.interp(i+1,j+1,k+1,shift)
                          - points.interp(i,  j+1,k,  shift)
                        )
                    );

                const vector aft =
                    0.5
                  * (
                        (
                            points.interp(i+1,j+1,k,shift)
                          - points.interp(i,  j,k,  shift)
                        )
                      ^ (
                            points.interp(i,  j+1,k,shift)
                          - points.interp(i+1,j,  k,shift)
                        )
                    );

                const vector fore =
                    0.5
                  * (
                        (
                            points.interp(i+1,j+1,k+1,shift)
                          - points.interp(i,  j,  k+1,shift)
                        )
                      ^ (
                            points.interp(i,  j+1,k+1,shift)
                          - points.interp(i+1,j,  k+1,shift)
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

    fn.correctBoundaryConditions();
    fa.correctBoundaryConditions();

    fan = fa*fn;

    fan.correctBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateCellCenters()
{
    meshField<vector,MeshType>& cc = cellCenters_;

    // First, set internal cell centers from the point coordinates

    forAll(cc, leveli)
    {
        const partLevelPoints& points = fvMsh_[leveli].points();

        forAll(cc[leveli], d)
        {
            meshDirection<vector,MeshType>& ccld = cc[leveli][d];

            const vector shift = MeshType::shift[d];

            forAllBlock(ccld, i, j, k)
            {
                ccld(i,j,k) =
                    0.125
                  * (
                        points.interp(i,  j,  k  ,shift)
                      + points.interp(i+1,j,  k  ,shift)
                      + points.interp(i,  j+1,k  ,shift)
                      + points.interp(i+1,j+1,k  ,shift)
                      + points.interp(i,  j,  k+1,shift)
                      + points.interp(i+1,j,  k+1,shift)
                      + points.interp(i,  j+1,k+1,shift)
                      + points.interp(i+1,j+1,k+1,shift)
                    );
            }
        }
    }

    // Next, project all inner cell centers along the point-to-point vector

    forAll(cc, leveli)
    {
        const partLevelPoints& points = fvMsh_[leveli].points();

        forAll(cc[leveli], d)
        {
            meshDirection<vector,MeshType>& ccld = cc[leveli][d];

            const vector shift = MeshType::shift[d];

            labelVector bo;

            for (bo.x() = -1; bo.x() <= 1; bo.x()++)
            for (bo.y() = -1; bo.y() <= 1; bo.y()++)
            for (bo.z() = -1; bo.z() <= 1; bo.z()++)
            if (cmptSum(cmptMag(bo)) > 0)
            {
                const labelVector S(ccld.boundaryStart(bo));
                const labelVector E(ccld.boundaryEnd(bo));

                const labelVector S2
                (
                    bo.x() == 0 ? 0 : (bo.x()+1)/2,
                    bo.y() == 0 ? 0 : (bo.y()+1)/2,
                    bo.z() == 0 ? 0 : (bo.z()+1)/2
                );

                const labelVector E2
                (
                    S2.x() + (bo.x() == 0 ? 2 : 1),
                    S2.y() + (bo.y() == 0 ? 2 : 1),
                    S2.z() + (bo.z() == 0 ? 2 : 1)
                );

                labelVector ijk, ijk2;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    vector p(0,0,0);
                    scalar c(0.0);

                    for (ijk2.x() = S2.x(); ijk2.x() < E2.x(); ijk2.x()++)
                    for (ijk2.y() = S2.y(); ijk2.y() < E2.y(); ijk2.y()++)
                    for (ijk2.z() = S2.z(); ijk2.z() < E2.z(); ijk2.z()++)
                    {
                        p += points.interp(ijk+ijk2, shift);
                        c += 1.0;
                    }

                    p /= c;

                    ccld(ijk+bo) = 2*p - ccld(ijk);
                }
            }
        }
    }

    // Don't call correctBoundaryConditions(), because this will update and
    // overwrite all boundaries. There's only a need to update parallel
    // boundaries.

    forAll(cc, l)
    {
        const label nReq = Pstream::nRequests();

        forAll(cc.boundaryConditions(), i)
        if (cc.boundaryConditions()[i].baseType() == PARALLELBC)
        {
            cc.boundaryConditions()[i].initEvaluate(l);
        }

        if (Pstream::parRun())
        {
            Pstream::waitRequests(nReq);
        }

        forAll(cc.boundaryConditions(), i)
        if (cc.boundaryConditions()[i].baseType() == PARALLELBC)
        {
            cc.boundaryConditions()[i].evaluate(l);
        }
    }
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateCellVolumes()
{
    meshField<scalar,MeshType>& cv = cellVolumes_;

    const meshField<faceVector,MeshType>& fn = faceNormals_;
    const meshField<faceScalar,MeshType>& fa = faceAreas_;
    const meshField<faceVector,MeshType>& fc = faceCenters_;

    forAll(cv, leveli)
    {
        // Cell volume is given by Wesseling's efficient formula (Wesseling, p.
        // 484)

        forAll(cv[leveli], d)
        {
            meshDirection<scalar,MeshType>& cvld = cv[leveli][d];

            const meshDirection<faceVector,MeshType>& fnld = fn[leveli][d];
            const meshDirection<faceScalar,MeshType>& fald = fa[leveli][d];
            const meshDirection<faceVector,MeshType>& fcld = fc[leveli][d];

            forAllBlock(cvld, i, j, k)
            {
                const vector Sx =
                    fnld(i,j,k).left() *fald(i,j,k).left()
                  + fnld(i,j,k).right()*fald(i,j,k).right();

                const vector Sy =
                    fnld(i,j,k).bottom()*fald(i,j,k).bottom()
                  + fnld(i,j,k).top()   *fald(i,j,k).top();

                const vector Sz =
                    fnld(i,j,k).aft() *fald(i,j,k).aft()
                  + fnld(i,j,k).fore()*fald(i,j,k).fore();

                const vector Dx = fcld(i,j,k).right() - fcld(i,j,k).left();
                const vector Dy = fcld(i,j,k).top()   - fcld(i,j,k).bottom();
                const vector Dz = fcld(i,j,k).fore()  - fcld(i,j,k).aft();

                cvld(i,j,k) = ((Dx & Sx) + (Dy & Sy) + (Dz & Sz))/6.0;
            }
        }
    }

    cv.correctBoundaryConditions();
}

template<class MeshType>
void fvMeshMetrics<MeshType>::calculateFaceDeltas()
{
    meshField<faceScalar,MeshType>& fd = faceDeltas_;
    const meshField<vector,MeshType>& cc = cellCenters_;

    forAll(cc, leveli)
    {
        forAll(cc[leveli], d)
        {
            meshDirection<faceScalar,MeshType>& fdld = fd[leveli][d];
            const meshDirection<vector,MeshType>& ccld = cc[leveli][d];

            forAllBlock(fdld, i, j, k)
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

    faceDeltas_.correctBoundaryConditions();
}

template<class MeshType>
fvMeshMetrics<MeshType>::fvMeshMetrics(const fvMesh& fvMsh)
:
    fvMsh_(fvMsh),
    faceCenters_
    (
        "faceCenters",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    faceNormals_
    (
        "faceNormals",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    faceAreas_
    (
        "faceCenters",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    faceAreaNormals_
    (
        "faceAreaNormals",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    cellCenters_
    (
        "cellCenters",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    cellVolumes_
    (
        "cellVolumes",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    faceDeltas_
    (
        "faceDeltas",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    )
{
    calculateFaceCenters();
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
