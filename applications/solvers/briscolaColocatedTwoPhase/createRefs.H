const colocatedFaceScalarField& fa =
    fvMsh.metrics<colocated>().faceAreas();

const colocatedLowerFaceScalarField& mu = twoPhase.mu<colocated>();

const colocatedScalarField& rho = twoPhase.rho<colocated>();
