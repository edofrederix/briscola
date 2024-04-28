#include "uniformMesh.H"
#include "rectilinearMesh.H"
#include "structuredMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

word directionText(labelVector v)
{
    v = cmptMag(v);

    if (cmptSum(v) == 0)
    {
        return word("no direction");
    }
    else if (cmptSum(v) == 1)
    {
        return
            word
            (
                v == unitX ? "x"
              : v == unitY ? "y"
              :              "z"
            )
          + word(" direction");
    }
    else if (cmptSum(v) == 2)
    {
        return
            word
            (
                v == unitXY ? "x and y"
              : v == unitXZ ? "x and z"
              :               "y and z"
            )
          + word(" directions");
    }
    else
    {
        return word("all directions");
    }
}

void generalInfo(const fvMesh& fvMsh)
{
    label nCells(cmptProduct(fvMsh[0].N()));
    reduce(nCells, sumOp<label>());

    Info<< "Global number of cells = " << nCells << endl;

    if (fvMsh.structured())
    {
        Info<< "Global mesh size = "
            << fvMsh.msh().cast<structuredMesh>().N() << endl;
    }

    Info<< "Global bounding box = " << fvMsh.msh().boundingBox() << endl;

    const colocatedScalarField& cv = fvMsh.metrics<colocated>().cellVolumes();

    Info<< "Global volume = " << gSum(cv)[0] << endl;

    Info<< "Global mesh is "
        << (fvMsh.structured() ? "structured" : "unstructured")
        << endl;

    if (fvMsh.structured())
    {
        Info<< "Global mesh is rectilinear in "
            << directionText(fvMsh.rectilinear()) << endl;
    }

    if (fvMsh.rectilinear() == unitXYZ)
    {
        const rectilinearMesh& rMesh = fvMsh.msh().cast<rectilinearMesh>();

        Info<< "    First base vector  = " << rMesh.base().x() << nl
            << "    Second base vector = " << rMesh.base().y() << nl
            << "    Third base vector  = " << rMesh.base().z() << endl;
    }

    if (fvMsh.structured())
    {
       Info<< "Global mesh is uniform in "
            << directionText(fvMsh.uniform()) << endl;
    }

    if (fvMsh.uniform() == unitXYZ)
    {
        const uniformMesh& uMesh = fvMsh.msh().cast<uniformMesh>();

        Info<< "    Uniform cell size = " << uMesh.cellSize() << endl;
    }

    Info<< "Brick topology is "
        << (fvMsh.msh().topology().structured() ? "structured" : "unstructured")
        << endl;

    if (fvMsh.structured())
    {
        Info<< "Brick topology is "
            << (fvMsh.msh().topology().aligned() ? "aligned" : "unaligned")
            << endl;

        Info<< "Brick topology is "
            << (fvMsh.msh().topology().rectilinear() ? "rectilinear" : "not rectilinear")
            << endl;
    }
}

void brickInfo(const fvMesh& fvMsh)
{
    Info<< endl;
    Info<< "Number of bricks = " << fvMsh.msh().bricks().size() << ":" << endl;

    forAll(fvMsh.msh().bricks(), bricki)
    {
        const brick& br = fvMsh.msh().bricks()[bricki];

        Info<< "    Brick " << br.num() << " has " << br.N() << " = "
            << cmptProduct(br.N())  << " cells" << endl
            << "        Vertex coordinates: " << endl
            << "            v0 = " << br.v0().v() << endl
            << "            v1 = " << br.v1().v() << endl
            << "            v2 = " << br.v2().v() << endl
            << "            v3 = " << br.v3().v() << endl
            << "            v4 = " << br.v4().v() << endl
            << "            v5 = " << br.v5().v() << endl
            << "            v6 = " << br.v6().v() << endl
            << "            v7 = " << br.v7().v() << endl;

        if (Pstream::parRun())
        {
            Info<< "        Decomposition vector = "
                << fvMsh.msh().decomp().decompPerBrick()[bricki] << endl;
        }
    }
}

void patchInfo(const fvMesh& fvMsh)
{
    Info<< endl;
    Info<< "Number of patches = "
        << fvMsh.msh().patches().size() << ":" << endl;

    forAll(fvMsh.msh().patches(), patchi)
    {
        const patch& p = fvMsh.msh().patches()[patchi];

        Info<< "    Patch " << p.name() << " is of type "
            << p.typeName() << " and has "
            << p.facePtrs().size() << " brick faces:" << endl;

        forAll(p.facePtrs(), facei)
        {
            const face& f = *p.facePtrs()[facei];

            Info<< "        Face " << f.num()
                << " of brick " << f.parentBrick().num()
                << endl;
        }
    }
}

void levelInfo(const fvMesh& fvMsh)
{
    Info<< endl;
    Info<< "Number of mesh levels = " << fvMsh.size() << ":" << endl;

    forAll(fvMsh, l)
    {
        labelVector N(fvMsh[l].N());

        List<labelVector> Ns(Pstream::nProcs(),N);

        Pstream::gatherList(Ns);

        Info<< "    Level " << l << ":" << endl;

        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            Info<< "        Processor " << proc << " has " << Ns[proc]
                << " = " << cmptProduct(Ns[proc]) << " cells" << endl;
        }
    }
}

void parallelInfo(const fvMesh& fvMsh)
{
    Info<< endl;
    Info<< "Number of processors = " << Pstream::nProcs() << ":" << endl;

    List<dictionary> pDicts(Pstream::nProcs(), dictionary());
    List<List<dictionary>> bDicts(Pstream::nProcs(), List<dictionary>(0));
    List<labelList> masterPatch(Pstream::nProcs(), labelList(0));

    const part& p = fvMsh[0];

    pDicts[Pstream::myProcNo()].add("rectilinear", p.rectilinear());
    pDicts[Pstream::myProcNo()].add("uniform", p.uniform());


    forAll(fvMsh.boundaries(), bi)
    {
        const boundary& b = fvMsh.boundaries()[bi];

        bDicts[Pstream::myProcNo()].append(b.dict());
        masterPatch[Pstream::myProcNo()].append(b.master());
    }

    Pstream::gatherList(pDicts);
    Pstream::gatherList(bDicts);
    Pstream::gatherList(masterPatch);

    if (Pstream::master())
    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        const label brickNum(fvMsh.msh().decomp().brickNumPerProc()[proc]);

        const labelVector brickPart
        (
            fvMsh.msh().decomp().brickPartPerProc()[proc]
        );

        const labelVector brickDecomp
        (
            fvMsh.msh().decomp().decompPerBrick()[brickNum]
        );

        const dictionary& pDict = pDicts[proc];

        const labelVector pSize
        (
            briscola::cmptDivide(fvMsh.msh().bricks()[brickNum].N(),brickDecomp)
        );

        Info<< "    Processor " << proc << " is in brick " << brickNum
            << " with part index = " << brickPart << endl
            << "        Number of cells = " << pSize << " = "
            << cmptProduct(pSize) << endl
            << "        Part is rectilinear in "
            << directionText(labelVector(pDict.lookup("rectilinear")))
            << endl
            << "        Part is uniform in "
            << directionText(labelVector(pDict.lookup("uniform")))
            << endl;

        Info<< "        Number of part boundaries = "
            << bDicts[proc].size();

        if (bDicts[proc].size() > 0)
        {
            Info<< ":" << endl;

            forAll(bDicts[proc], bi)
            {
                const dictionary& dict = bDicts[proc][bi];

                const word name(dict.lookup("name"));
                const word type(dict.lookup("type"));
                const labelVector bo(dict.lookup("offset"));
                const labelTensor T(dict.lookup("T"));

                Info<< "                Boundary " << name << nl
                    << "                type   = " << type << nl
                    << "                bo     = " << bo << nl
                    << "                T      = " << T << nl
                    << "                master = " << masterPatch[proc][bi] << endl;
            }
        }
        else
        {
            Info<< endl;
        }
    }
}

void parallelConnectivityInfo(const fvMesh& fvMsh)
{
    Info<< endl;

    List<List<dictionary>> bDicts(Pstream::nProcs());

    forAll(fvMsh.boundaries(), bi)
    {
        const boundary& b = fvMsh.boundaries()[bi];

        if (b.castable<parallelBoundary>())
        {
            dictionary dict(b.dict());

            dict.add("master", b.master());
            dict.add("tag", b.cast<parallelBoundary>().tag());

            bDicts[Pstream::myProcNo()].append
            (
                dictionary(dict)
            );
        }
    }

    Pstream::gatherList(bDicts);

    label count = 0;

    forAll(bDicts, proci)
    {
        count += bDicts[proci].size();
    }

    Info<< "Found " << count << " boundaries that require communication,"
        << " distributed across " << Pstream::nProcs()
        << " processors";

    if (Pstream::master())
    {
        Info<< ":" << endl;

        forAll(bDicts, proci)
        forAll(bDicts[proci], bi)
        {
            const dictionary& dict0 = bDicts[proci][bi];

            const labelVector offset0
            (
                dict0.lookup("offset")
            );

            const labelVector neighborOffset0
            (
                dict0.lookup("neighborOffset")
            );

            const label procj =
                readLabel(dict0.lookup("neighborProcNum"));

            const word name0(dict0.lookup("name"));
            const labelTensor T0(dict0.lookup("T"));
            const label tag0(readLabel(dict0.lookup("tag")));
            const label master0(readLabel(dict0.lookup("master")));

            // Find the opposite of this parallel boundary

            bool found = false;

            forAll(bDicts[procj], bj)
            {
                const dictionary& dict1 = bDicts[procj][bj];

                const labelVector offset1
                (
                    dict1.lookup("offset")
                );

                const labelVector neighborOffset1
                (
                    dict1.lookup("neighborOffset")
                );

                const word name1(dict1.lookup("name"));
                const labelTensor T1(dict1.lookup("T"));
                const label tag1(readLabel(dict1.lookup("tag")));
                const label master1(readLabel(dict1.lookup("master")));

                if
                (
                    offset0 == neighborOffset1
                 && neighborOffset0 == offset1
                )
                {
                    found = true;

                    Info<< "    boundary " << name0 << " of processor " << proci
                        << " is connected with boundary " << name1
                        << " of processor " << procj << endl;

                    Info<< "        at processor " << proci << nl
                        << "            bo     = " << offset0 << nl
                        << "            T      = " << T0 << nl
                        << "            tag    = " << tag0 << nl
                        << "            master = " << master0 << endl;

                    Info<< "        at processor " << procj << nl
                        << "            bo     = " << offset1 << nl
                        << "            T      = " << T1 << nl
                        << "            tag    = " << tag1 << nl
                        << "            master = " << master1 << endl;

                    if (master0 == master1)
                    {
                        FatalErrorInFunction
                            << "Mismatch in master/slave. Both boundaries are "
                            << ((master0 == 0) ? "slave" : "master")
                            << endl << exit(FatalError);
                    }

                    break;
                }
            }

            if (!found)
            {
                FatalErrorInFunction
                    << "Could not find matching processor boundary for "
                    << nl << dict0 << nl << " of processor " << proci
                    << " in processor " << procj
                    << ". Available boundaries are:" << nl;

                forAll(bDicts[procj], bj)
                if
                (
                    word(bDicts[procj][bj].lookup("type"))
                == "parallel"
                )
                {
                    FatalErrorInFunction << bDicts[procj][bj] << nl;
                }

                FatalErrorInFunction << endl << exit(FatalError);
            }
        }
    }
    else
    {
        Info<< endl;
    }
}

}

}

}