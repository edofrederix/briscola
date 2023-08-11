IOdictionary vofDict
(
    IOobject
    (
        "briscolaVofDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

vof vf(vofDict, fvMsh);
