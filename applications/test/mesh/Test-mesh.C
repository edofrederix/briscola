#include "arguments.H"
#include "Time.H"
#include "uniformMesh.H"
#include "lineEdge.H"

using namespace Foam;
using namespace briscola;

void checkBrick(const brick& b)
{
    // Check faces

    for (int i = 0; i < 6; i++)
    {
        const labelVector offset = faceOffsets[i];

        labelBlock v(b.v());

        for (int d = 0; d < 3; d++)
            if (offset[d] != 0)
                v = v.slice((offset[d]+1)/2, d);

        const face& f = b.faces()[i];

        for (int j = 0; j < 4; j++)
            if (f.v()(j) != v(j))
                FatalErrorInFunction
                    << "Test 1a failed" << endl << abort(FatalError);

        if (b.f(i) != f)
            FatalErrorInFunction
                << "Test 1b failed" << endl << abort(FatalError);
    }

    // Check edges

    for (int i = 0; i < 12; i++)
    {
        const labelVector offset = edgeOffsets[i];

        labelBlock v(b.v());

        for (int d = 0; d < 3; d++)
            if (offset[d] != 0)
                v = v.slice((offset[d]+1)/2, d);

        const edge& e = b.e(i);

        if
        (
           !(
                (e.v()(0) == v(0) && e.v()(1) == v(1))
             || (e.v()(1) == v(0) && e.v()(0) == v(1))
            )
        )
        {
            FatalErrorInFunction
                << "Test 1c failed" << endl << abort(FatalError);
        }
    }

    // Check vertices

    for (int i = 0; i < 8; i++)
    {
        const labelVector offset = vertexOffsets[i];

        labelBlock v(b.v());

        for (int d = 0; d < 3; d++)
            if (offset[d] != 0)
                v = v.slice((offset[d]+1)/2, d);

        const vertex& vert = b.v(i);

        if (vert.j() != v(0))
        {
            FatalErrorInFunction
                << "Test 1d failed" << endl << abort(FatalError);
        }
    }

    // Check volume

    if (b.volume() != (b.rightHanded() ? 1.0 : -1.0) && Pstream::master())
    {
        Info << b.volume() << " " << endl;
        FatalErrorInFunction
           << "Test 1e failed" << endl << abort(FatalError);
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

    autoPtr<mesh> mshPtr(mesh::New(meshDict));

    // Check upcasts. Should all work because mesh is uniform.

    if (mshPtr->castable<uniformMesh>())
    {
        mshPtr->cast<uniformMesh>();
    }
    else
    {
        FatalErrorInFunction
            << "Test 1a failed" << endl << abort(FatalError);
    }

    if (mshPtr->castable<rectilinearMesh>())
    {
        mshPtr->cast<rectilinearMesh>();
    }
    else
    {
        FatalErrorInFunction
            << "Test 1b failed" << endl << abort(FatalError);
    }

    if (mshPtr->castable<structuredMesh>())
    {
        mshPtr->cast<structuredMesh>();
    }
    else
    {
        FatalErrorInFunction
            << "Test 1c failed" << endl << abort(FatalError);
    }

    if (mshPtr->castable<unstructuredMesh>())
    {
        mshPtr->cast<unstructuredMesh>();
    }
    else
    {
        FatalErrorInFunction
            << "Test 1d failed" << endl << abort(FatalError);
    }

    const uniformMesh& msh = mshPtr->cast<uniformMesh>();

    // Mesh should be structured, rectilinear in three directions and uniform in
    // three directions

    if (!msh.structured())
    {
        FatalErrorInFunction
            << "Test 2 failed" << endl << abort(FatalError);
    }

    if (msh.rectilinear() != unitXYZ)
    {
        FatalErrorInFunction
            << "Test 3 failed" << endl << abort(FatalError);
    }

    if (msh.uniform() != unitXYZ)
    {
        FatalErrorInFunction
            << "Test 4 failed" << endl << abort(FatalError);
    }

    // Copy constructors

    brick b(msh.bricks()[0]);
    face f(b.f0());
    lineEdge e(f.e0());
    vertex v(e.v0());

    if (b.f0() != f)
        FatalErrorInFunction
            << "Test 5 failed" << endl << abort(FatalError);

    if (f.e0() != e)
        FatalErrorInFunction
            << "Test 6 failed" << endl << abort(FatalError);

    if (e.v0() != v)
        FatalErrorInFunction
            << "Test 7 failed" << endl << abort(FatalError);

    // Check bricks and all possible transformations of bricks

    forAll(msh.bricks(), i)
    {
        const brick& b0 = msh.bricks()[i];

        checkBrick(b0);

        // Check rotations (includes eye rotation)

        for (int j = 0; j < 12; j++)
        {
            brick b(b0);
            b.transform(rotations[j/4][j%4]);
            checkBrick(b);

            // Check reflections of rotations

            for (int j = 0; j < 3; j++)
            {
                brick b2(b);
                b2.transform(reflections[j]);
                checkBrick(b2);
            }

            // Check permutations of rotations

            for (int j = 0; j < 3; j++)
            {
                brick b2(b);
                b2.transform(permutations[j]);
                checkBrick(b2);
            }

            // Check permutations of reflections of rotations

            for (int j = 0; j < 3; j++)
            {
                brick b2(b);
                b2.transform(reflections[j]);
                checkBrick(b2);

                for (int k = 0; k < 3; k++)
                {
                    b2.transform(permutations[k]);
                    checkBrick(b2);
                }
            }

            // Check reflections of permutations of rotations

            for (int j = 0; j < 3; j++)
            {
                brick b2(b);
                b2.transform(permutations[j]);
                checkBrick(b2);

                for (int k = 0; k < 3; k++)
                {
                    b2.transform(reflections[k]);
                    checkBrick(b2);
                }
            }
        }

        // Check rotate member function

        for (int d = 0; d < 3; d++)
        {
            for (int j = 0; j < 4; j++)
            {
                brick b1(b0);
                brick b2(b0);

                const labelTensor T = rotations[d][j];

                b1.transform(T);
                b2.rotate(j,d);

                checkBrick(b1);
                checkBrick(b2);

                if (b1 != b2)
                    FatalErrorInFunction
                        << "Test 8 failed" << endl << abort(FatalError);
            }
        }
    }

    // Brick topology should be structured, aligned and rectilinear

    const brickTopology& topo = msh.topology();

    if (!topo.structured())
    {
        FatalErrorInFunction
            << "Test 9 failed" << endl << abort(FatalError);
    }

    if (!topo.aligned())
    {
        FatalErrorInFunction
            << "Test 10 failed" << endl << abort(FatalError);
    }

    if (!topo.rectilinear())
    {
        FatalErrorInFunction
            << "Test 11 failed" << endl << abort(FatalError);
    }

    // Because the brick topology is aligned, brick links should have identity
    // transformations

    forAll(topo.links(), bricki)
    {
        forAll(topo.links()[bricki].faceLinks(), linki)
        if (topo.links()[bricki].faceLinks().set(linki))
        {
            const brickFaceLink& link =
                topo.links()[bricki].faceLinks()[linki];

            if (link.T() != eye)
                FatalErrorInFunction
                    << "Test 12 failed" << endl << abort(FatalError);
        }

        forAll(topo.links()[bricki].edgeLinks(), linki)
        if (topo.links()[bricki].edgeLinks().set(linki))
        {
            const brickEdgeLink& link =
                topo.links()[bricki].edgeLinks()[linki];

            if (link.T() != eye)
                FatalErrorInFunction
                    << "Test 13 failed" << endl << abort(FatalError);
        }

        forAll(topo.links()[bricki].vertexLinks(), linki)
        if (topo.links()[bricki].vertexLinks().set(linki))
        {
            const brickVertexLink& link =
                topo.links()[bricki].vertexLinks()[linki];

            if (link.T() != eye)
                FatalErrorInFunction
                    << "Test 14 failed" << endl << abort(FatalError);
        }
    }

    const scalarList& dx = msh.cellSizes()[0];
    const scalarList& dy = msh.cellSizes()[1];
    const scalarList& dz = msh.cellSizes()[2];

    forAll(dx, i)
    {
        if (Foam::mag(dx[i]-1.0/24) > 1e-12)
            FatalErrorInFunction
                << "Test 15a failed" << endl << abort(FatalError);
    }

    forAll(dy, i)
    {
        if (Foam::mag(dy[i]-1.0/16) > 1e-12)
            FatalErrorInFunction
                << "Test 15b failed" << endl << abort(FatalError);
    }

    forAll(dz, i)
    {
        if (Foam::mag(dz[i]-1.0/16) > 1e-12)
            FatalErrorInFunction
                << "Test 15c failed" << endl << abort(FatalError);
    }
}
