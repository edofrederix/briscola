#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"
#include "IO.H"

#include "restrictionSchemes.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

const vector L(1, 1.5, 4);

template<class MeshType>
void testCellCenters(const fvMesh& fvMsh)
{
    const meshField<vector,MeshType>& c =
        fvMsh.metrics<MeshType>().cellCenters();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.decomp().myBrickDecomp()))
    );

    forAll(c, l)
    forAll(c[l], d)
    for (label i = c[l][d].S().x()-1; i < c[l][d].E().x()+1; i++)
    for (label j = c[l][d].S().y()-1; j < c[l][d].E().y()+1; j++)
    for (label k = c[l][d].S().z()-1; k < c[l][d].E().z()+1; k++)
    {
        const vector cc
        (
            cmptMultiply
            (
                cmptDivide
                (
                    vector(i+0.5, j+0.5, k+0.5)
                  + MeshType::shift[d],
                    vector(fvMsh[l].N())
                )
              + vector(fvMsh.decomp().myBrickPart()),
                Lp
            )
        );

        if (mag(c[l][d](i,j,k) - cc) > 1e-12)
        {
            FatalErrorInFunction
                << "test 1 failed" << abort(FatalError);
        }
    }
}

template<class MeshType>
void testCellVolumes(const fvMesh& fvMsh)
{
    const meshField<scalar,MeshType>& v =
        fvMsh.metrics<MeshType>().cellVolumes();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.decomp().myBrickDecomp()))
    );

    forAll(v, l)
    {
        const scalar V
        (
            Lp.x()/fvMsh[l].N().x()
          * Lp.y()/fvMsh[l].N().y()
          * Lp.z()/fvMsh[l].N().z()
        );

        forAll(v[l], d)
        forAllCells(v[l][d], i, j, k)
        {
            if (mag(v[l][d](i,j,k) - V) > 1e-12)
            {
                FatalErrorInFunction
                    << "test 2 failed" << abort(FatalError);
            }
        }
    }
}

template<class MeshType>
void testFaceCenters(const fvMesh& fvMsh)
{
    const meshField<faceVector,MeshType>& c =
        fvMsh.metrics<MeshType>().faceCenters();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.decomp().myBrickDecomp()))
    );

    forAll(c, l)
    forAll(c[l], d)
    forAllCells(c[l][d], i, j, k)
    for (label o = 0; o < 6; o++)
    {
        const vector cc
        (
            cmptMultiply
            (
                cmptDivide
                (
                    vector(i+0.5, j+0.5, k+0.5)
                  + MeshType::shift[d]
                  + vector(faceOffsets[o])*0.5,
                    vector(fvMsh[l].N())
                )
                + vector(fvMsh.decomp().myBrickPart()),
                Lp
            )
        );

        if (mag(c[l][d](i,j,k)[o] - cc) > 1e-12)
        {
            FatalErrorInFunction
                << "test 3 failed" << abort(FatalError);
        }
    }
}

template<class MeshType>
void testFaceAreas(const fvMesh& fvMsh)
{
    const meshField<faceScalar,MeshType>& a =
        fvMsh.metrics<MeshType>().faceAreas();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.decomp().myBrickDecomp()))
    );

    forAll(a, l)
    {
        const vector A
        (
            Lp.y()/fvMsh[l].N().y() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.y()/fvMsh[l].N().y()
        );

        forAll(a[l], d)
        forAllCells(a[l][d], i, j, k)
        for (label o = 0; o < 6; o++)
        {
            if (mag(a[l][d](i,j,k)[o] - A[o/2]) > 1e-12)
            {
                FatalErrorInFunction
                    << "test 4 failed" << abort(FatalError);
            }
        }
    }
}

template<class MeshType>
void testFaceNormals(const fvMesh& fvMsh)
{
    const meshField<faceVector,MeshType>& n =
        fvMsh.metrics<MeshType>().faceNormals();

    forAll(n, l)
    {
        forAll(n[l], d)
        forAllCells(n[l][d], i, j, k)
        for (label o = 0; o < 6; o++)
        {
            if (mag(n[l][d](i,j,k)[o] - vector(units[o/2])) > 1e-12)
            {
                FatalErrorInFunction
                    << "test 5 failed" << abort(FatalError);
            }
        }
    }
}

template<class MeshType>
void testFaceAreaNormals(const fvMesh& fvMsh)
{
    const meshField<faceVector,MeshType>& n =
        fvMsh.metrics<MeshType>().faceAreaNormals();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.decomp().myBrickDecomp()))
    );

    forAll(n, l)
    {
        const vector A
        (
            Lp.y()/fvMsh[l].N().y() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.y()/fvMsh[l].N().y()
        );

        forAll(n[l], d)
        forAllCells(n[l][d], i, j, k)
        for (label o = 0; o < 6; o++)
        {
            if (mag(n[l][d](i,j,k)[o] - A[o/2]*vector(units[o/2])) > 1e-12)
            {
                FatalErrorInFunction
                    << "test 6 failed" << abort(FatalError);
            }
        }
    }
}

template<class MeshType>
void testFaceDeltas(const fvMesh& fvMsh)
{
    const meshField<faceScalar,MeshType>& fd =
        fvMsh.metrics<MeshType>().faceDeltas();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.decomp().myBrickDecomp()))
    );

    forAll(fd, l)
    {
        const vector D
        (
            1.0/(Lp.x()/fvMsh[l].N().x()),
            1.0/(Lp.y()/fvMsh[l].N().y()),
            1.0/(Lp.z()/fvMsh[l].N().z())
        );

        forAll(fd[l], d)
        forAllCells(fd[l][d], i, j, k)
        for (label o = 0; o < 6; o++)
        {
            if (mag(fd[l][d](i,j,k)[o] - D[o/2]) > 1e-12)
            {
                FatalErrorInFunction
                    << "test 7 failed" << abort(FatalError);
            }
        }
    }
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"

    IOdictionary meshDict
    (
        IOobject
        (
            runTime.system()/"briscolaMeshDict",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    fvMesh fvMsh(meshDict, runTime);

    #include "createBriscolaIO.H"

    // Test cell centers

    testCellCenters<colocated>(fvMsh);
    testCellCenters<staggered>(fvMsh);

    // Test cell volumes

    testCellVolumes<colocated>(fvMsh);
    testCellVolumes<staggered>(fvMsh);

    // Test face centers

    testFaceCenters<colocated>(fvMsh);
    testFaceCenters<staggered>(fvMsh);

    // Test face areas

    testFaceAreas<colocated>(fvMsh);
    testFaceAreas<staggered>(fvMsh);

    // Test face normals

    testFaceNormals<colocated>(fvMsh);
    testFaceNormals<staggered>(fvMsh);

    // Test face area normals

    testFaceAreaNormals<colocated>(fvMsh);
    testFaceAreaNormals<staggered>(fvMsh);

    // Test face deltas

    testFaceDeltas<colocated>(fvMsh);
    testFaceDeltas<staggered>(fvMsh);
}
