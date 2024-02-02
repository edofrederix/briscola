const colocatedFaceScalarField& fa =
    fvMsh.metrics<colocated>().faceAreas();

const colocatedLowerFaceScalarField& mu = tpm.mu<colocated>();
