namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
void checkParallelFaceCenters(const fvMesh& fvMsh)
{
    const meshField<faceVector,MeshType>& fc =
        fvMsh.metrics<MeshType>().faceCenters();

    for (label i = 0; i < 6; i++)
    {
        const labelVector bo(faceOffsets[i]);

        bool found = false;

        forAll(fc.boundaryConditions(), bci)
        {
            const boundaryCondition<faceVector,MeshType>& bc =
                fc.boundaryConditions()[bci];

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
    const meshField<faceScalar,MeshType>& fd =
        fvMsh.metrics<MeshType>().faceDeltas();

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
                        << "Found mismatch in face centers at "
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