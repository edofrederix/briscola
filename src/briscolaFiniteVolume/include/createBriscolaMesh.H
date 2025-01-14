IOdictionary meshDict
(
    IOobject
    (
        "briscolaMeshDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

fvMesh fvMsh(meshDict, runTime);
