const colocatedFaceScalarField& fa =
    fvMsh.metrics<colocated>().faceAreas();

const staggeredScalarField& rho = tpm.rho<staggered>();
const colocatedScalarField& coloRho = tpm.rho<colocated>();
const staggeredFaceScalarField& mu = tpm.mu<staggered>();
const staggeredScalarField minRho = tpm.minRho<staggered>().ref();
const colocatedScalarField coloMinRho = tpm.minRho<colocated>().ref();
const staggeredScalarField meanRho = tpm.meanRho<staggered>().ref();

const staggeredScalarField minRhoInv("minRhoInv", 1.0/minRho);
const colocatedFaceScalarField coloMinRhof("coloMinRhof",ex::interp(coloMinRho));