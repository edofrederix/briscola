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
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAll(c, l)
    forAll(c[l], d)
    for (label i = c.I(l,d).left()  -1; i < c.I(l,d).right()+1; i++)
    for (label j = c.I(l,d).bottom()-1; j < c.I(l,d).top()  +1; j++)
    for (label k = c.I(l,d).aft()   -1; k < c.I(l,d).fore() +1; k++)
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
              + vector(fvMsh.msh().decomp().myBrickPart()),
                Lp
            )
        );

        if (mag(c(l,d,i,j,k) - cc) > 1e-12)
            FatalErrorInFunction
                << "test 1 failed" << abort(FatalError);
    }
}

template<class MeshType>
void testCellVolumes(const fvMesh& fvMsh)
{
    const meshField<scalar,MeshType>& v =
        fvMsh.metrics<MeshType>().cellVolumes();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllCells(v, l, d, i, j, k)
    {
        const scalar V
        (
            Lp.x()/fvMsh[l].N().x()
          * Lp.y()/fvMsh[l].N().y()
          * Lp.z()/fvMsh[l].N().z()
        );

        if (mag(v(l,d,i,j,k) - V) > 1e-12)
            FatalErrorInFunction
                << "test 2 failed" << abort(FatalError);
    }
}

template<class MeshType>
void testFaceCenters(const fvMesh& fvMsh)
{
    const meshField<faceVector,MeshType>& c =
        fvMsh.metrics<MeshType>().faceCenters();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllCells(c, l, d, i, j, k)
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
              + vector(fvMsh.msh().decomp().myBrickPart()),
                Lp
            )
        );

        if (mag(c(l,d,i,j,k)[o] - cc) > 1e-12)
            FatalErrorInFunction
                << "test 3a failed" << abort(FatalError);
    }
}

template<class MeshType>
void testEdgeCenters(const fvMesh& fvMsh)
{
    const meshField<edgeVector,MeshType>& c =
        fvMsh.metrics<MeshType>().edgeCenters();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllCells(c, l, d, i, j, k)
    for (label o = 0; o < 12; o++)
    {
        const vector cc
        (
            cmptMultiply
            (
                cmptDivide
                (
                    vector(i+0.5, j+0.5, k+0.5)
                  + MeshType::shift[d]
                  + vector(edgeOffsets[o])*0.5,
                    vector(fvMsh[l].N())
                )
              + vector(fvMsh.msh().decomp().myBrickPart()),
                Lp
            )
        );

        if (mag(c(l,d,i,j,k)[o] - cc) > 1e-12)
            FatalErrorInFunction
                << "test 3b failed" << abort(FatalError);
    }
}

template<class MeshType>
void testVertexCenters(const fvMesh& fvMsh)
{
    const meshField<vertexVector,MeshType>& c =
        fvMsh.metrics<MeshType>().vertexCenters();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllCells(c, l, d, i, j, k)
    for (label o = 0; o < 8; o++)
    {
        const vector cc
        (
            cmptMultiply
            (
                cmptDivide
                (
                    vector(i+0.5, j+0.5, k+0.5)
                  + MeshType::shift[d]
                  + vector(vertexOffsets[o])*0.5,
                    vector(fvMsh[l].N())
                )
              + vector(fvMsh.msh().decomp().myBrickPart()),
                Lp
            )
        );

        if (mag(c(l,d,i,j,k)[o] - cc) > 1e-12)
            FatalErrorInFunction
                << "test 3c failed" << abort(FatalError);
    }
}

template<class MeshType>
void testFaceAreas(const fvMesh& fvMsh)
{
    const meshField<faceScalar,MeshType>& a =
        fvMsh.metrics<MeshType>().faceAreas();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllCells(a, l, d, i, j, k)
    {
        const vector A
        (
            Lp.y()/fvMsh[l].N().y() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.y()/fvMsh[l].N().y()
        );

        for (label o = 0; o < 6; o++)
            if (mag(a(l,d,i,j,k)[o] - A[o/2]) > 1e-12)
                FatalErrorInFunction
                    << "test 4 failed" << abort(FatalError);
    }
}

template<class MeshType>
void testFaceNormals(const fvMesh& fvMsh)
{
    const meshField<faceVector,MeshType>& fn =
        fvMsh.metrics<MeshType>().faceNormals();

    forAllCells(fn, l, d, i, j, k)
    for (label o = 0; o < 6; o++)
    {
        label lr = 2*(o%2)-1;

        if (mag(fn(l,d,i,j,k)[o] - lr*vector(briscola::units[o/2])) > 1e-12)
            FatalErrorInFunction
                << "test 5 failed" << abort(FatalError);
    }

    // Check if face normals are pointing outwards

    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    const meshField<faceVector,MeshType>& fc =
        fvMsh.metrics<MeshType>().faceCenters();

    forAllCells(fn, l, d, i, j, k)
    for (label o = 0; o < 6; o++)
    {
        vector f = fc(l,d,i,j,k)[o] - cc(l,d,i,j,k);
        vector n = fn(l,d,i,j,k)[o];

        if ((f & n) < 0)
            FatalErrorInFunction
                << "test 6 failed" << abort(FatalError);
    }
}

template<class MeshType>
void testFaceAreaNormals(const fvMesh& fvMsh)
{
    const meshField<faceVector,MeshType>& fan =
        fvMsh.metrics<MeshType>().faceAreaNormals();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllCells(fan, l, d, i, j, k)
    {
        const vector A
        (
            Lp.y()/fvMsh[l].N().y() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.y()/fvMsh[l].N().y()
        );

        for (label o = 0; o < 6; o++)
        {
            label lr = 2*(o%2)-1;

            if
            (
                mag(fan(l,d,i,j,k)[o] - lr*A[o/2]*vector(briscola::units[o/2]))
              > 1e-12
            )
            {
                FatalErrorInFunction
                    << "test 7 failed" << abort(FatalError);
            }
        }
    }

    // Check if face area normals are pointing outwards

    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    const meshField<faceVector,MeshType>& fc =
        fvMsh.metrics<MeshType>().faceCenters();

    forAllCells(fan, l, d, i, j, k)
    for (label o = 0; o < 6; o++)
    {
        vector f = fc(l,d,i,j,k)[o] - cc(l,d,i,j,k);
        vector n = fan(l,d,i,j,k)[o];

        if ((f & n) < 0)
            FatalErrorInFunction
                << "test 8 failed" << abort(FatalError);
    }
}

template<class MeshType>
void testFaceDeltas(const fvMesh& fvMsh)
{
    const meshField<faceScalar,MeshType>& fd =
        fvMsh.metrics<MeshType>().faceDeltas();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllCells(fd, l, d, i, j, k)
    {
        const vector D
        (
            1.0/(Lp.x()/fvMsh[l].N().x()),
            1.0/(Lp.y()/fvMsh[l].N().y()),
            1.0/(Lp.z()/fvMsh[l].N().z())
        );

        for (label o = 0; o < 6; o++)
            if (mag(fd(l,d,i,j,k)[o] - D[o/2]) > 1e-12)
                FatalErrorInFunction
                    << "test 9 failed" << abort(FatalError);
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

    // Test edge centers

    testEdgeCenters<colocated>(fvMsh);
    testEdgeCenters<staggered>(fvMsh);

    // Test vertex centers

    testVertexCenters<colocated>(fvMsh);
    testVertexCenters<staggered>(fvMsh);

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
