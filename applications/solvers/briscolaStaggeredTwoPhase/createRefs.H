const colocatedFaceScalarField& fa =
    fvMsh.metrics<colocated>().faceAreas();

const staggeredScalarField& rho = tpm.rho<staggered>();
const staggeredLowerFaceScalarField& mu = tpm.mu<staggered>();
