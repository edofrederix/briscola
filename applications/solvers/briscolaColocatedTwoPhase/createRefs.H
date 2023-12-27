const colocatedFaceScalarField& fa =
    fvMsh.metrics<colocated>().faceAreas();

const colocatedScalarField& rho = tpm.rho<colocated>();
const colocatedLowerFaceScalarField& mu = tpm.mu<colocated>();
