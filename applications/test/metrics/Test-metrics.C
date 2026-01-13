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

        if (mag(c(l,d,i,j,k) - cc)/mag(L) > 1e-12)
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

        // Slightly larger tolerance because of larger roundoff errors in the
        // volume computation
        if (mag(v(l,d,i,j,k) - V)/V > 1e-11)
            FatalErrorInFunction
                << "test 2 failed" << abort(FatalError);
    }
}

template<class MeshType>
void testFaceCenters(const fvMesh& fvMsh)
{
    const faceField<vector,MeshType>& fc =
        fvMsh.metrics<MeshType>().faceCenters();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllFaces(fc, fd, l, d, i, j, k)
    {
        const vector cc
        (
            cmptMultiply
            (
                cmptDivide
                (
                    vector(i+0.5, j+0.5, k+0.5)
                  + MeshType::shift[d]
                  + vector(faceOffsets[fd*2])*0.5,
                    vector(fvMsh[l].N())
                )
              + vector(fvMsh.msh().decomp().myBrickPart()),
                Lp
            )
        );

        if (mag(fc[fd](l,d,i,j,k) - cc)/mag(L) > 1e-12)
            FatalErrorInFunction
                << "test 3a failed" << abort(FatalError);
    }
}

template<class MeshType>
void testEdgeCenters(const fvMesh& fvMsh)
{
    const meshField<edgeVector,MeshType> ec
    (
        fvMsh.metrics<MeshType>().aos().edgeCenters()
    );

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllCells(ec, l, d, i, j, k)
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

        if (mag(ec(l,d,i,j,k)[o] - cc)/mag(L) > 1e-12)
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

        if (mag(c(l,d,i,j,k)[o] - cc)/mag(L) > 1e-12)
            FatalErrorInFunction
                << "test 3c failed" << abort(FatalError);
    }
}

template<class MeshType>
void testFaceAreas(const fvMesh& fvMsh)
{
    const faceField<scalar,MeshType>& fa =
        fvMsh.metrics<MeshType>().faceAreas();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllFaces(fa, fd, l, d, i, j, k)
    {
        const vector A
        (
            Lp.y()/fvMsh[l].N().y() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.y()/fvMsh[l].N().y()
        );

        if (mag(fa[fd](l,d,i,j,k) - A[fd])/A[fd] > 1e-12)
            FatalErrorInFunction
                << "test 4 failed" << abort(FatalError);
    }
}

template<class MeshType>
void testFaceNormals(const fvMesh& fvMsh)
{
    const faceField<vector,MeshType>& fn =
        fvMsh.metrics<MeshType>().faceNormals();

    forAllFaces(fn, fd, l, d, i, j, k)
    {
        if (mag(fn[fd](l,d,i,j,k) + vector(briscola::units[fd])) > 1e-12)
            FatalErrorInFunction
                << "test 5 failed" << abort(FatalError);
    }

    // Check if face normals are pointing outwards

    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    const faceField<vector,MeshType>& fc =
        fvMsh.metrics<MeshType>().faceCenters();

    forAllFaces(fn, fd, l, d, i, j, k)
    {
        vector f = fc[fd](l,d,i,j,k) - cc(l,d,i,j,k);
        vector n = fn[fd](l,d,i,j,k);

        if ((f & n) < 0)
            FatalErrorInFunction
                << "test 6 failed" << abort(FatalError);
    }
}

template<class MeshType>
void testFaceAreaNormals(const fvMesh& fvMsh)
{
    const faceField<vector,MeshType>& fan =
        fvMsh.metrics<MeshType>().faceAreaNormals();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllFaces(fan, fd, l, d, i, j, k)
    {
        const vector A
        (
            Lp.y()/fvMsh[l].N().y() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.z()/fvMsh[l].N().z(),
            Lp.x()/fvMsh[l].N().x() * Lp.y()/fvMsh[l].N().y()
        );

        if
        (
            mag(fan[fd](l,d,i,j,k) + A[fd]*vector(briscola::units[fd]))
          / A[fd]
          > 1e-12
        )
        {
            FatalErrorInFunction
                << "test 7 failed" << abort(FatalError);
        }
    }

    // Check if face area normals are pointing outwards

    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    const faceField<vector,MeshType>& fc =
        fvMsh.metrics<MeshType>().faceCenters();

    forAllFaces(fan, fd, l, d, i, j, k)
    {
        vector f = fc[fd](l,d,i,j,k) - cc(l,d,i,j,k);
        vector n = fan[fd](l,d,i,j,k);

        if ((f & n) < 0)
            FatalErrorInFunction
                << "test 8 failed" << abort(FatalError);
    }
}

template<class MeshType>
void testFaceDeltas(const fvMesh& fvMsh)
{
    const faceField<scalar,MeshType>& delta =
        fvMsh.metrics<MeshType>().faceDeltas();

    const vector Lp
    (
        cmptDivide(L, vector(fvMsh.msh().decomp().myBrickDecomp()))
    );

    forAllFaces(delta, fd, l, d, i, j, k)
    {
        const vector D
        (
            1.0/(Lp.x()/fvMsh[l].N().x()),
            1.0/(Lp.y()/fvMsh[l].N().y()),
            1.0/(Lp.z()/fvMsh[l].N().z())
        );

        if (mag(delta[fd](l,d,i,j,k) - D[fd])/D[fd] > 1e-12)
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
