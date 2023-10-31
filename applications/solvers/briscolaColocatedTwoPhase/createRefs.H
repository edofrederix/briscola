const colocatedFaceScalarField& fa =
    fvMsh.metrics<colocated>().faceAreas();

const colocatedScalarField& rho = tpm.rho<colocated>();
const colocatedFaceScalarField& mu = tpm.mu<colocated>();