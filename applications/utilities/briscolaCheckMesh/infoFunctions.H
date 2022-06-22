namespace Foam
{

namespace briscola
{

namespace fv
{

void generalInfo(const fvMesh& fvMsh)
{
    label nCells(cmptProduct(fvMsh.N()));
    reduce(nCells, sumOp<label>());

    Info<< "Number of cells = " << nCells << endl;

    const colocatedScalarField& cv = fvMsh.metrics<colocated>().cellVolumes();
    Info<< "Mesh volume = " << gSum(cv)[0] << endl;

    Info<< "Bricks are "
        << (fvMsh.topology().structured() ? "structured" : "unstructured")
        << endl;

    Info<< endl;
}

void brickInfo(const fvMesh& fvMsh)
{
    Info<< "Number of bricks = " << fvMsh.bricks().size() << ":" << endl;

    forAll(fvMsh.bricks(), bricki)
    {
        const brick& br = fvMsh.bricks()[bricki];

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
            Info<< "        Decomposition vector = " << fvMsh.decomp().decompPerBrick()[bricki] << endl;
        }
    }

    Info<< endl;
}

void patchInfo(const fvMesh& fvMsh)
{
    Info<< "Number of patches = " << fvMsh.patches().size() << ":" << endl;

    forAll(fvMsh.patches(), patchi)
    {
        const patch& p = fvMsh.patches()[patchi];

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

    Info<< endl;
}

void parallelInfo(const fvMesh& fvMsh)
{
    Info<< "Number of processors = " << Pstream::nProcs() << ":" << endl;

    List<List<dictionary>> dicts(Pstream::nProcs(), List<dictionary>(0));

    forAll(fvMsh.partPatches(), patchi)
    {
        const partPatch& patch = fvMsh.partPatches()[patchi];

        dicts[Pstream::myProcNo()].append(patch.dict());
    }

    Pstream::gatherList(dicts);

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        const label brickNum(fvMsh.decomp().brickNumPerProc()[proc]);
        const labelVector brickPart(fvMsh.decomp().brickPartPerProc()[proc]);
        const labelVector brickDecomp(fvMsh.decomp().decompPerBrick()[brickNum]);

        Info<< "    Processor " << proc << " is in brick " << brickNum
            << " with part index = " << brickPart << endl
            << "        Number of cells = "
            << briscola::cmptDivide(fvMsh.bricks()[brickNum].N(),brickDecomp) << " = "
            << cmptProduct(briscola::cmptDivide(fvMsh.bricks()[brickNum].N(),brickDecomp)) << endl;

        Info<< "        Number of part patches = "
            << dicts[proc].size();

        if (dicts[proc].size() > 0)
        {
            Info<< ":" << endl;

            forAll(dicts[proc], patchi)
            {
                const dictionary& dict = dicts[proc][patchi];

                const word name(dict.lookup("name"));
                const word type(dict.lookup("type"));
                const labelVector bo(dict.lookup("boundaryOffset"));
                const labelTensor T(dict.lookup("T"));

                Info<< "            Part patch " << name << nl
                    << "                type   = " << type << nl
                    << "                bo     = " << bo << nl
                    << "                offset = " << bo << nl
                    << "                T      = " << T << endl;
            }
        }
        else
        {
            Info<< endl;
        }
    }

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

    Info<< endl;
}

void parallelConnectivityInfo(const fvMesh& fvMsh)
{
    List<List<dictionary>> patchDicts(Pstream::nProcs());

    forAll(fvMsh.partPatches(), patchi)
    {
        const partPatch& patch = fvMsh.partPatches()[patchi];

        if
        (
            patch.type() == "parallel"
        ||  patch.type() == "periodic"
        )
        {
            patchDicts[Pstream::myProcNo()].append
            (
                dictionary(patch.dict())
            );
        }
    }

    Pstream::gatherList(patchDicts);

    label count = 0;

    forAll(patchDicts, proci)
    {
        count += patchDicts[proci].size();
    }

    Info<< "Found " << count << " patches that require communication,"
        << " distributed across " << Pstream::nProcs()
        << " processors:" << endl;

    if (Pstream::master())
    {
        forAll(patchDicts, proci)
        forAll(patchDicts[proci], patchi)
        {
            const dictionary& dict0 = patchDicts[proci][patchi];

            const labelVector boundaryOffset0
            (
                dict0.lookup("boundaryOffset")
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

            // Find the opposite of this parallel patch

            bool found = false;

            forAll(patchDicts[procj], patchj)
            {
                const dictionary& dict1 = patchDicts[procj][patchj];

                const labelVector boundaryOffset1
                (
                    dict1.lookup("boundaryOffset")
                );

                const labelVector neighborOffset1
                (
                    dict1.lookup("neighborOffset")
                );

                const word name1(dict1.lookup("name"));
                const labelTensor T1(dict1.lookup("T"));
                const label tag1(readLabel(dict1.lookup("tag")));

                if
                (
                    boundaryOffset0 == neighborOffset1
                 && neighborOffset0 == boundaryOffset1
                )
                {
                    found = true;

                    Info<< "    patch " << name0 << " of processor " << proci
                        << " is connected with patch " << name1 << " of processor "
                        << procj << endl;

                    Info<< "        at processor " << proci << nl
                        << "            offset = " << boundaryOffset0 << nl
                        << "            T      = " << T0 << nl
                        << "            tag    = " << tag0 << endl;

                    Info<< "        at processor " << procj << nl
                        << "            offset = " << boundaryOffset1 << nl
                        << "            T      = " << T1 << nl
                        << "            tag    = " << tag1 << endl;

                    break;
                }
            }

            if (!found)
            {
                FatalErrorInFunction
                    << "Could not find matching processor patch for patch "
                    << nl << dict0 << nl << " of processor " << proci
                    << " in processor " << procj << ". Available patches are:"
                    << nl;

                forAll(patchDicts[procj], patchj)
                if
                (
                    word(patchDicts[procj][patchj].lookup("type"))
                == "parallel"
                )
                {
                    FatalErrorInFunction << patchDicts[procj][patchj] << nl;
                }

                FatalErrorInFunction << endl << exit(FatalError);
            }
        }

        Info<< endl;
    }
}

}

}

}