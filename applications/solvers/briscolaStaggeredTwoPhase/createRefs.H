const staggeredScalarField& rho = tpm.rho<staggered>();
const colocatedScalarField& coloRho = tpm.rho<colocated>();
const staggeredFaceScalarField& mu = tpm.mu<staggered>();
const staggeredScalarField minRho = tpm.minRho<staggered>().ref();

const staggeredScalarField minRhoInv("minRhoInv", 1.0/minRho);
const staggeredFaceScalarField minRhoInvf("minRhoInv",1.0/ex::interp(minRho));
const staggeredFaceScalarField minRhof("minRhof",ex::interp(minRho));