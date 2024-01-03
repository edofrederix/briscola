namespace Foam
{

namespace briscola
{

namespace fv
{

void checkParallelPatches(const fvMesh& fvMsh)
{
    const mesh& msh = fvMsh.msh();
    const decomposition& decomp = msh.decomp();

    const faceLabel& faceNeighbors =
        decomp.faceNeighborsPerProc()[Pstream::myProcNo()];

    const edgeLabel& edgeNeighbors =
        decomp.edgeNeighborsPerProc()[Pstream::myProcNo()];

    const vertexLabel& vertexNeighbors =
        decomp.vertexNeighborsPerProc()[Pstream::myProcNo()];

    forAll(msh.boundaries(), i)
    {
        const boundary& b = msh.boundaries()[i];

        if (b.castable<parallelBoundary>())
        {
            const labelVector bo = b.offset();

            if (b.offsetDegree() == 1)
            {
                const label faceNum = faceNumber(bo);

                if (faceNeighbors[faceNum] != b.neighborProcNum())
                    FatalErrorInFunction
                        << "Mismatch in face boundary neighbor "
                        << "processor number."
                        << endl << exit(FatalError);
            }
            else if (b.offsetDegree() == 2)
            {
                const label edgeNum = edgeNumber(bo);

                if (edgeNeighbors[edgeNum] != b.neighborProcNum())
                    FatalErrorInFunction
                        << "Mismatch in edge boundary neighbor "
                        << "processor number."
                        << endl << exit(FatalError);
            }
            else
            {
                const label vertexNum = vertexNumber(bo);

                if (vertexNeighbors[vertexNum] != b.neighborProcNum())
                    FatalErrorInFunction
                        << "Mismatch in vertex boundary neighbor "
                        << "processor number."
                        << endl << exit(FatalError);
            }
        }
    }
}

template<class MeshType>
void checkInternalFaceCenters(const fvMesh& fvMsh)
{
    const meshField<faceVector,MeshType>& fc =
        fvMsh.metrics<MeshType>().faceCenters();

    forAllCells(fc, l, d, i, j, k)
    {
        const labelVector N(fc.N(l,d));
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
                        fc(l,d,ijk0)[face0]
                      - fc(l,d,ijk1)[face1]
                    )
                  > 1e-12
                )
                {
                    FatalErrorInFunction
                        << "Found mismatch in face centers at "
                        << "offset = " << offset << nl
                        << "mesh type = " << MeshType::typeName << nl
                        << "mesh direction = " << d << nl
                        << "first cell face = " << fc(l,d,ijk0)[face0] << nl
                        << "first cell faces = " << fc(l,d,ijk0) << nl
                        << "second cell face = "
                        << fc(l,d,ijk1)[face1] << nl
                        << "second cell faces = " << fc(l,d,ijk1)<< nl
                        << "level = " << l << nl
                        << "first cell index = " << ijk0 << nl
                        << "second cell index = " << ijk1 << nl
                        << "direction size = " << fc.N(l,d) << endl
                        << exit(FatalError);
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

    forAllCells(fd, l, d, i, j, k)
    {
        const labelVector N(fd.N(l,d));

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
                        1.0/fd(l,d,ijk0)[face0]
                      - 1.0/fd(l,d,ijk1)[face1]
                    )
                  > 1e-12
                )
                {
                    FatalErrorInFunction
                        << "Found mismatch in face deltas at "
                        << "offset = " << offset << nl
                        << "mesh type = " << MeshType::typeName << nl
                        << "mesh direction = " << d << nl
                        << "first cell delta = "
                        << fd(l,d,ijk0)[face0] << nl
                        << "first cell deltas = " << fd(l,d,ijk0) << nl
                        << "second cell delta = "
                        << fd(l,d,ijk1)[face1] << nl
                        << "second cell deltas = " << fd(l,d,ijk1)<< nl
                        << "level = " << l << nl
                        << "fist cell index = " << ijk0 << nl
                        << "second cell index = " << ijk1 << nl
                        << "direction size = " << fd.N(l,d) << endl
                        << exit(FatalError);
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

            if (bc.offset() == bo && bc.baseType() == PARALLELBC)
            {
                found = true;
                slave = bc.mshBoundary().slave();
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
                    Info
                        << "Found mismatch in face centers at "
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
                bc.offset() == bo
             && (
                    bc.baseType() == PERIODICBC
                 || bc.baseType() == PARALLELBC
                )
            )
            {
                found = true;
                slave = bc.mshBoundary().slave();
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
                        << "Found mismatch in face deltas at "
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
