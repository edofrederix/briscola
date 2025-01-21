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

autoPtr<TwoPhaseModel<colocated>> twoPhasePtr
(
    TwoPhaseModel<colocated>::New(fvMsh, twoPhaseDict)
);

TwoPhaseModel<colocated>& twoPhase = twoPhasePtr();

const Switch reduced = twoPhase.reduced();
