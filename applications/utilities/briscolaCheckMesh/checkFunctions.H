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

        forAllBlock(fcld, i, j, k)
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
                            << "fc0[face] = " << fcld(ijk0)[face0] << nl
                            << "fc0 = " << fcld(ijk0) << nl
                            << "fc1[face] = " << fcld(ijk1)[face1] << nl
                            << "fc1 = " << fcld(ijk1)<< nl
                            << "level = " << l << nl
                            << "ijk0 = " << ijk0 << endl
                            << "ijk1 = " << ijk1 << endl
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

        forAllBlock(fdld, i, j, k)
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
                            << "Found mismatch in face delas at "
                            << "offset = " << offset << nl
                            << "mesh type = " << MeshType::typeName << nl
                            << "mesh direction = " << d << nl
                            << "fd0[face] = " << fdld(ijk0)[face0] << nl
                            << "fd0 = " << fdld(ijk0) << nl
                            << "fd1[face] = " << fdld(ijk1)[face1] << nl
                            << "fd1 = " << fdld(ijk1)<< nl
                            << "level = " << l << nl
                            << "ijk0 = " << ijk0 << endl
                            << "ijk1 = " << ijk1 << endl
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

        forAll(fc.boundaryConditions(), bci)
        {
            const boundaryCondition<faceVector,MeshType>& bc =
                fc.boundaryConditions()[bci];

            if (bc.boundaryOffset() == bo && bc.baseType() == PARALLELBC)
            {
                found = true;
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

            const labelVector S(fcld.boundaryStart(bo));
            const labelVector E(fcld.boundaryEnd(bo));

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
                    FatalErrorInFunction
                        << "Found mismatch in face centers at "
                        << "boundary offset = " << bo << nl
                        << "mesh type = " << MeshType::typeName << nl
                        << "mesh direction = " << d << nl
                        << "fc1[face] = " << fcld(ijk)[dir1] << nl
                        << "fc1 = " << fcld(ijk) << nl
                        << "fc2[face] = " << fcld(ijk+bo)[dir2] << nl
                        << "fc2 = " << fcld(ijk+bo)<< nl
                        << "level = " << l << nl
                        << "ijk = " << ijk << endl
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

            const labelVector S(fdld.boundaryStart(bo));
            const labelVector E(fdld.boundaryEnd(bo));

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
                        << "fd1[face] = " << fdld(ijk)[dir1] << nl
                        << "fd1 = " << fdld(ijk) << nl
                        << "fd2[face] = " << fdld(ijk+bo)[dir2] << nl
                        << "fd2 = " << fdld(ijk+bo)<< nl
                        << "level = " << l << nl
                        << "ijk = " << ijk << endl
                        << exit(FatalError);
                }
            }
        }
    }
}

}

}

}