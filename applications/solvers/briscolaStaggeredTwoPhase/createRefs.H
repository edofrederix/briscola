const colocatedFaceScalarField& fa =
    fvMsh.metrics<colocated>().faceAreas();

const staggeredScalarField& rho = tpm.rho<staggered>();
const colocatedScalarField& coloRho = tpm.rho<colocated>();
const staggeredFaceScalarField& mu = tpm.mu<staggered>();