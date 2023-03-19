namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
void checkInternalFaceCenters(const fvMesh& fvMsh)
{
    const meshField<faceVector,MeshType>& fc =
        fvMsh.metrics<MeshType>().faceCenters();

    forAll(fc, l)
    forAll(fc[l], d)
    {
        const meshDirection<faceVector,MeshType>& fcld = fc[l][d];

        const labelVector N(fcld.N());

        forAllCells(fcld, i, j, k)
        {
            const labelVector ijk0(i,j,k);

            for (label face0 = 0; face0 < 6; face0++)
            {
                const labelVector offset(faceOffsets[face0]);
                const labelVector ijk1(ijk0+offset);

                const label face1 = faceNumber(-offset);

                if
                (
                    ijk1.x() >= 0 && ijk1.x() < N.x()
                 && ijk1.y() >= 0 && ijk1.y() < N.y()
                 && ijk1.z() >= 0 && ijk1.z() < N.z()
                )
                {
                    if
                    (
                        Foam::mag
                        (
                            fcld(ijk0)[face0]
                          - fcld(ijk1)[face1]
                        )
                      > 1e-12
                    )
                    {
                        FatalErrorInFunction
                            << "Found mismatch in face centers at "
                            << "offset = " << offset << nl
                            << "mesh type = " << MeshType::typeName << nl
                            << "mesh direction = " << d << nl
                            << "first cell face = " << fcld(ijk0)[face0] << nl
                            << "first cell faces = " << fcld(ijk0) << nl
                            << "second cell face = " << fcld(ijk1)[face1] << nl
                            << "second cell faces = " << fcld(ijk1)<< nl
                            << "level = " << l << nl
                            << "first cell index = " << ijk0 << nl
                            << "second cell index = " << ijk1 << nl
                            << "direction size = " << fcld.N() << endl
                            << exit(FatalError);
                    }
                }
            }
        }

    }
}

template<class MeshType>
void checkInternalFaceDeltas(const fvMesh& fvMsh)
{
    const meshField<faceScalar,MeshType>& fd =
        fvMsh.metrics<MeshType>().faceDeltas();

    forAll(fd, l)
    forAll(fd[l], d)
    {
        const meshDirection<faceScalar,MeshType>& fdld = fd[l][d];

        const labelVector N(fdld.N());

        forAllCells(fdld, i, j, k)
        {
            const labelVector ijk0(i,j,k);

            for (label face0 = 0; face0 < 6; face0++)
            {
                const labelVector offset(faceOffsets[face0]);
                const labelVector ijk1(ijk0+offset);

                const label face1 = faceNumber(-offset);

                if
                (
                    ijk1.x() >= 0 && ijk1.x() < N.x()
                 && ijk1.y() >= 0 && ijk1.y() < N.y()
                 && ijk1.z() >= 0 && ijk1.z() < N.z()
                )
                {
                    if
                    (
                        Foam::mag
                        (
                            1.0/fdld(ijk0)[face0]
                          - 1.0/fdld(ijk1)[face1]
                        )
                      > 1e-12
                    )
                    {
                        FatalErrorInFunction
                            << "Found mismatch in face deltas at "
                            << "offset = " << offset << nl
                            << "mesh type = " << MeshType::typeName << nl
                            << "mesh direction = " << d << nl
                            << "first cell delta = " << fdld(ijk0)[face0] << nl
                            << "first cell deltas = " << fdld(ijk0) << nl
                            << "second cell delta = " << fdld(ijk1)[face1] << nl
                            << "second cell deltas = " << fdld(ijk1)<< nl
                            << "level = " << l << nl
                            << "fist cell index = " << ijk0 << nl
                            << "second cell index = " << ijk1 << nl
                            << "direction size = " << fdld.N() << endl
                            << exit(FatalError);
                    }
                }
            }
        }

    }
}

template<class MeshType>
void checkParallelFaceCenters(const fvMesh& fvMsh)
{
    // We need to correct boundary conditions for this check, because by default
    // ghost cell face centers are zero

    meshField<faceVector,MeshType>& fc =
        const_cast<meshField<faceVector,MeshType>&>
        (
            fvMsh.metrics<MeshType>().faceCenters()
        );

    fc.correctBoundaryConditions();

    for (label i = 0; i < 6; i++)
    {
        const labelVector bo(faceOffsets[i]);

        bool found = false;
        bool slave = false;

        forAll(fc.boundaryConditions(), bci)
        {
            const boundaryCondition<faceVector,MeshType>& bc =
                fc.boundaryConditions()[bci];

            if (bc.boundaryOffset() == bo && bc.baseType() == PARALLELBC)
            {
                found = true;
                slave = bc.slave();
                break;
            }
        }

        if (!found)
            continue;

        const label dir1(faceNumber(bo));
        const label dir2(faceNumber(-bo));

        forAll(fc, l)
        forAll(fc[l], d)
        {
            const meshDirection<faceVector,MeshType>& fcld = fc[l][d];

            labelVector S(fcld.activeBoundaryStart(bo));
            labelVector E(fcld.activeBoundaryEnd(bo));

            labelVector ijk;

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                if
                (
                    Foam::mag
                    (
                        fcld(ijk+bo)[dir2]
                      - fcld(ijk)[dir1]
                    )
                  > 1e-12
                )
                {
                    Info
                        << "Found mismatch in face centers at "
                        << "boundary offset = " << bo << nl
                        << "mesh type = " << MeshType::typeName << nl
                        << "mesh direction = " << d << nl
                        << "internal cell face = " << fcld(ijk)[dir1] << nl
                        << "internal cell faces = " << fcld(ijk) << nl
                        << "ghost cell face = " << fcld(ijk+bo)[dir2] << nl
                        << "ghost cell face vector = " << fcld(ijk+bo)<< nl
                        << "level = " << l << nl
                        << "internal cell index = " << ijk << nl
                        << "ghost cell index = " << ijk+bo << nl
                        << "direction size = " << fcld.N() << nl
                        << "shifted = " << Switch(fcld.shifted(bo)) << nl
                        << "slave = " << Switch(slave) << endl
                        << exit(FatalError);
                }
            }
        }
    }
}

template<class MeshType>
void checkParallelFaceDeltas(const fvMesh& fvMsh)
{
    // We need to correct boundary conditions for this check, because by default
    // ghost cell face deltas are zero

    meshField<faceScalar,MeshType>& fd =
        const_cast<meshField<faceScalar,MeshType>&>
        (
            fvMsh.metrics<MeshType>().faceDeltas()
        );

    fd.correctBoundaryConditions();

    for (label i = 0; i < 6; i++)
    {
        const labelVector bo(faceOffsets[i]);

        bool found = false;
        bool slave = false;

        forAll(fd.boundaryConditions(), bci)
        {
            const boundaryCondition<faceScalar,MeshType>& bc =
                fd.boundaryConditions()[bci];

            if
            (
                bc.boundaryOffset() == bo
             && (
                    bc.baseType() == PERIODICBC
                 || bc.baseType() == PARALLELBC
                )
            )
            {
                found = true;
                slave = bc.slave();
                break;
            }
        }

        if (!found)
            continue;

        const label dir1(faceNumber(bo));
        const label dir2(faceNumber(-bo));

        forAll(fd, l)
        forAll(fd[l], d)
        {
            const meshDirection<faceScalar,MeshType>& fdld = fd[l][d];

            labelVector S(fdld.activeBoundaryStart(bo));
            labelVector E(fdld.activeBoundaryEnd(bo));

            labelVector ijk;

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                if
                (
                    Foam::mag
                    (
                        1.0/fdld(ijk+bo)[dir2]
                      - 1.0/fdld(ijk)[dir1]
                    )
                    > 1e-12
                )
                {
                    FatalErrorInFunction
                        << "Found mismatch in face deltas at "
                        << "boundary offset = " << bo << nl
                        << "mesh direction = " << d << nl
                        << "internal cell delta = " << fdld(ijk)[dir1] << nl
                        << "internal cell deltas = " << fdld(ijk) << nl
                        << "ghost cell delta = " << fdld(ijk+bo)[dir2] << nl
                        << "ghost cell deltas = " << fdld(ijk+bo)<< nl
                        << "level = " << l << nl
                        << "internal cell index = " << ijk << nl
                        << "ghost cell index = " << ijk+bo << nl
                        << "direction size = " << fdld.N() << endl
                        << "shifted = " << Switch(fdld.shifted(bo)) << nl
                        << "slave = " << Switch(slave) << endl
                        << exit(FatalError);
                }
            }
        }
    }
}

}

}

}