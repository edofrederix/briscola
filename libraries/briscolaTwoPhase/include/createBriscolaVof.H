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

autoPtr<vof> vfPtr(vof::New(vofDict, fvMsh));

vof& vf = vfPtr();
