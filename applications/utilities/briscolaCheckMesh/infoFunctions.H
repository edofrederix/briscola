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

    if (fvMsh.topology().structured())
    {
        Info<< "Global mesh shape = " << fvMsh.N() << endl;
    }

    const colocatedScalarField& cv = fvMsh.metrics<colocated>().cellVolumes();

    Info<< "Global volume = " << gSum(cv)[0] << endl;

    Info<< "Global mesh is "
        << (fvMsh.structured() ? "structured" : "unstructured")
        << endl;

    if (fvMsh.structured())
    {
        Info<< "Global mesh is rectilinear in "
            << directionText(fvMsh.rectilinear()) << endl;

        Info<< "Global mesh is uniform in "
            << directionText(fvMsh.uniform()) << endl;
    }

    Info<< "Brick topology is "
        << (fvMsh.topology().structured() ? "structured" : "unstructured")
        << endl;

    if (fvMsh.topology().structured())
    {
        Info<< "Brick topology is "
            << (fvMsh.topology().aligned() ? "aligned" : "unaligned")
            << endl;

        Info<< "Brick topology is "
            << (fvMsh.topology().rectilinear() ? "rectilinear" : "not rectilinear")
            << endl;
    }
}

void brickInfo(const fvMesh& fvMsh)
{
    Info<< endl;
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
}

void patchInfo(const fvMesh& fvMsh)
{
    Info<< endl;
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

    List<dictionary> partDicts(Pstream::nProcs(), dictionary());
    List<List<dictionary>> patchDicts(Pstream::nProcs(), List<dictionary>(0));
    List<labelList> masterPatch(Pstream::nProcs(), labelList(0));

    const partLevel& part = fvMsh[0];

    partDicts[Pstream::myProcNo()].add("rectilinear", part.rectilinear());
    partDicts[Pstream::myProcNo()].add("uniform", part.uniform());


    forAll(fvMsh.partPatches(), patchi)
    {
        const partPatch& patch = fvMsh.partPatches()[patchi];

        patchDicts[Pstream::myProcNo()].append(patch.dict());
        masterPatch[Pstream::myProcNo()].append(patch.master());
    }

    Pstream::gatherList(partDicts);
    Pstream::gatherList(patchDicts);
    Pstream::gatherList(masterPatch);

    if (Pstream::master())
    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        const label brickNum(fvMsh.decomp().brickNumPerProc()[proc]);
        const labelVector brickPart(fvMsh.decomp().brickPartPerProc()[proc]);
        const labelVector brickDecomp(fvMsh.decomp().decompPerBrick()[brickNum]);

        const dictionary& partDict = partDicts[proc];

        const labelVector partSize
        (
            briscola::cmptDivide(fvMsh.bricks()[brickNum].N(),brickDecomp)
        );

        Info<< "    Processor " << proc << " is in brick " << brickNum
            << " with part index = " << brickPart << endl
            << "        Number of cells = " << partSize << " = "
            << cmptProduct(partSize) << endl
            << "        Part is rectilinear in "
            << directionText(labelVector(partDict.lookup("rectilinear")))
            << endl
            << "        Part is uniform in "
            << directionText(labelVector(partDict.lookup("uniform")))
            << endl;

        Info<< "        Number of part patches = "
            << patchDicts[proc].size();

        if (patchDicts[proc].size() > 0)
        {
            Info<< ":" << endl;

            forAll(patchDicts[proc], patchi)
            {
                const dictionary& dict = patchDicts[proc][patchi];

                const word name(dict.lookup("name"));
                const word type(dict.lookup("type"));
                const labelVector bo(dict.lookup("boundaryOffset"));
                const labelTensor T(dict.lookup("T"));

                Info<< "            Part patch " << name << nl
                    << "                type   = " << type << nl
                    << "                bo     = " << bo << nl
                    << "                offset = " << bo << nl
                    << "                T      = " << T << nl
                    << "                master = " << masterPatch[proc][patchi] << endl;
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
            dictionary dict(patch.dict());

            dict.add("master", patch.master());

            patchDicts[Pstream::myProcNo()].append
            (
                dictionary(dict)
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
        << " processors";

    if (Pstream::master())
    {
        Info<< ":" << endl;

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
            const label master0(readLabel(dict0.lookup("master")));

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
                const label master1(readLabel(dict1.lookup("master")));

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
                        << "            tag    = " << tag0 << nl
                        << "            master = " << master0 << endl;

                    Info<< "        at processor " << procj << nl
                        << "            offset = " << boundaryOffset1 << nl
                        << "            T      = " << T1 << nl
                        << "            tag    = " << tag1 << nl
                        << "            master = " << master1 << endl;

                    if (master0 == master1)
                    {
                        FatalErrorInFunction
                            << "Mismatch in master/slave. Both patches are "
                            << ((master0 == 0) ? "slave" : "master")
                            << endl << exit(FatalError);
                    }

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
    }
    else
    {
        Info<< endl;
    }
}

}

}

}