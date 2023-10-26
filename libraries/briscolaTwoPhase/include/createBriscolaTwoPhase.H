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

autoPtr<twoPhaseModel> tpmPtr(twoPhaseModel::New(twoPhaseDict, fvMsh));

twoPhaseModel& tpm = tpmPtr();
