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

autoPtr<twoPhaseModel> twoPhasePtr(twoPhaseModel::New(fvMsh, twoPhaseDict));

twoPhaseModel& twoPhase = twoPhasePtr();
