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

    meshField<faceVector,MeshType> f
    (
        fvMsh.metrics<MeshType>().aos().faceCenters()
    );

    f.correctBoundaryConditions();

    for (label i = 0; i < 6; i++)
    {
        const labelVector bo(faceOffsets[i]);

        const label d1(faceNumber(bo));
        const label d2(faceNumber(-bo));

        forAll(f, l)
        {
            bool found = false;
            bool slave = false;

            forAll(f[l].boundaryConditions(), bci)
            {
                const boundaryCondition<faceVector,MeshType>& bc =
                    f[l].boundaryConditions()[bci];

                if (bc.offset() == bo && bc.baseType() == PARALLELBC)
                {
                    found = true;
                    slave = bc.b().template cast<parallelBoundary>().slave();

                    break;
                }
            }

            if (!found)
                continue;

            forAll(f[l], d)
            {
                const meshDirection<faceVector,MeshType>& fd = f[l][d];

                labelVector S(f.fvMsh().template S<MeshType>(l,d,bo));
                labelVector E(f.fvMsh().template E<MeshType>(l,d,bo));

                labelVector ijk;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                    if(Foam::mag(fd(ijk+bo)[d2] - fd(ijk)[d1]) > 1e-12)
                        FatalErrorInFunction
                            << "Found mismatch in face centers at parallel "
                            << "boundary offset = " << bo << nl
                            << "mesh type = " << MeshType::typeName << nl
                            << "mesh direction = " << d << nl
                            << "internal cell face = " << fd(ijk)[d1] << nl
                            << "internal cell faces = " << fd(ijk) << nl
                            << "ghost cell face = " << fd(ijk+bo)[d2] << nl
                            << "ghost cell face vector = " << fd(ijk+bo)<< nl
                            << "level = " << l << nl
                            << "internal cell index = " << ijk << nl
                            << "ghost cell index = " << ijk+bo << nl
                            << "direction size = " << f.N(l,d) << nl
                            << "shifted = " << Switch(MeshType::shifted(d,bo))
                            << nl << "slave = " << Switch(slave) << endl
                            << abort(FatalError);
            }
        }
    }
}

template<class MeshType>
void checkParallelFaceDeltas(const fvMesh& fvMsh)
{
    // We need to correct boundary conditions for this check, because by default
    // ghost cell face deltas are zero

    meshField<faceScalar,MeshType> f
    (
        fvMsh.metrics<MeshType>().aos().faceDeltas()
    );

    f.correctBoundaryConditions();

    for (label i = 0; i < 6; i++)
    {
        const labelVector bo(faceOffsets[i]);

        const label d1(faceNumber(bo));
        const label d2(faceNumber(-bo));

        forAll(f, l)
        {
            bool found = false;
            bool slave = false;

            forAll(f[l].boundaryConditions(), bci)
            {
                const boundaryCondition<faceScalar,MeshType>& bc =
                    f[l].boundaryConditions()[bci];

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

            forAll(f[l], d)
            {
                const meshDirection<faceScalar,MeshType>& fd = f[l][d];

                labelVector S(f.fvMsh().template S<MeshType>(l,d,bo));
                labelVector E(f.fvMsh().template E<MeshType>(l,d,bo));

                labelVector ijk;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                    if (Foam::mag(1.0/fd(ijk+bo)[d2]-1.0/fd(ijk)[d1]) > 1e-12)
                        FatalErrorInFunction
                            << "Found mismatch in face deltas at parallel "
                            << "boundary offset = " << bo << nl
                            << "mesh direction = " << d << nl
                            << "internal cell delta = " << fd(ijk)[d1] << nl
                            << "internal cell deltas = " << fd(ijk) << nl
                            << "ghost cell delta = " << fd(ijk+bo)[d2] << nl
                            << "ghost cell deltas = " << fd(ijk+bo)<< nl
                            << "level = " << l << nl
                            << "internal cell index = " << ijk << nl
                            << "ghost cell index = " << ijk+bo << nl
                            << "direction size = " << f.N(l,d) << endl
                            << "shifted = " << Switch(MeshType::shifted(d,bo))
                            << nl << "slave = " << Switch(slave) << endl
                            << abort(FatalError);
            }
        }
    }
}

}

}

}
