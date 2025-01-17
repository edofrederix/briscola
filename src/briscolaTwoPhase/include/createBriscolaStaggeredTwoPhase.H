IOdictionary twoPhaseDict
(
    IOobject
    (
        "briscolaTwoPhaseDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<TwoPhaseModel<staggered>> twoPhasePtr
(
    TwoPhaseModel<staggered>::New(fvMsh, twoPhaseDict)
);

TwoPhaseModel<staggered>& twoPhase = twoPhasePtr();

const Switch reduced = twoPhase.reduced();
