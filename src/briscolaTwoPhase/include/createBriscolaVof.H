colocatedScalarField alpha
(
    "alpha",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true,
    false
);

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

autoPtr<normalScheme> normalSchemePtr
(
    normalScheme::New(fvMsh, vofDict.subDict("normalScheme"), alpha).ptr()
);

normalScheme& normal = normalSchemePtr();

autoPtr<vof> vfPtr(vof::New(fvMsh, vofDict, normal, alpha));

vof& vf = vfPtr();
