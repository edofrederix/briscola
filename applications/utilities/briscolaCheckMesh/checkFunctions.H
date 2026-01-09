namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
void checkParallelFaceCenters(const fvMesh& fvMsh)
{
    // We need to correct boundary conditions for this check, because by default
    // ghost cell face centers are zero

    meshField<faceVector,MeshType> fc
    (
        fvMsh.metrics<MeshType>().aos().faceCenters()
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

            if (bc.offset() == bo && bc.baseType() == PARALLELBC)
            {
                found = true;
                slave = bc.b().template cast<parallelBoundary>().slave();

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
            labelVector S(fc.fvMsh().template S<MeshType>(l,d,bo));
            labelVector E(fc.fvMsh().template E<MeshType>(l,d,bo));

            labelVector ijk;

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                if
                (
                    Foam::mag
                    (
                        fc(l,d,ijk+bo)[dir2]
                      - fc(l,d,ijk)[dir1]
                    )
                  > 1e-12
                )
                {
                    FatalErrorInFunction
                        << "Found mismatch in face centers at parallel "
                        << "boundary offset = " << bo << nl
                        << "mesh type = " << MeshType::typeName << nl
                        << "mesh direction = " << d << nl
                        << "internal cell face = " << fc(l,d,ijk)[dir1] << nl
                        << "internal cell faces = " << fc(l,d,ijk) << nl
                        << "ghost cell face = " << fc(l,d,ijk+bo)[dir2] << nl
                        << "ghost cell face vector = " << fc(l,d,ijk+bo)<< nl
                        << "level = " << l << nl
                        << "internal cell index = " << ijk << nl
                        << "ghost cell index = " << ijk+bo << nl
                        << "direction size = " << fc.N(l,d) << nl
                        << "shifted = " << Switch(MeshType::shifted(d,bo))
                        << nl
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

    meshField<faceScalar,MeshType> fd
    (
        fvMsh.metrics<MeshType>().aos().faceDeltas()
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
                bc.offset() == bo
             && (
                    bc.baseType() == PERIODICBC
                 || bc.baseType() == PARALLELBC
                )
            )
            {
                found = true;
                slave = bc.b().template cast<parallelBoundary>().slave();

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
            labelVector S(fd.fvMsh().template S<MeshType>(l,d,bo));
            labelVector E(fd.fvMsh().template E<MeshType>(l,d,bo));

            labelVector ijk;

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                if
                (
                    Foam::mag
                    (
                        1.0/fd(l,d,ijk+bo)[dir2]
                      - 1.0/fd(l,d,ijk)[dir1]
                    )
                    > 1e-12
                )
                {
                    FatalErrorInFunction
                        << "Found mismatch in face deltas at parallel "
                        << "boundary offset = " << bo << nl
                        << "mesh direction = " << d << nl
                        << "internal cell delta = " << fd(l,d,ijk)[dir1] << nl
                        << "internal cell deltas = " << fd(l,d,ijk) << nl
                        << "ghost cell delta = " << fd(l,d,ijk+bo)[dir2] << nl
                        << "ghost cell deltas = " << fd(l,d,ijk+bo)<< nl
                        << "level = " << l << nl
                        << "internal cell index = " << ijk << nl
                        << "ghost cell index = " << ijk+bo << nl
                        << "direction size = " << fd.N(l,d) << endl
                        << "shifted = " << Switch(MeshType::shifted(d,bo))
                        << nl
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
