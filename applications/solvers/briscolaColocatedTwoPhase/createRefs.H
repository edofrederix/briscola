const colocatedFaceScalarField& fa =
    fvMsh.metrics<colocated>().faceAreas();

const colocatedScalarField& rho = tpm.rho<colocated>();
const colocatedLowerFaceScalarField& mu = tpm.mu<colocated>();
const colocatedScalarField minRho = tpm.minRho<colocated>().ref();

const colocatedScalarField minRhoInv("minRhoInv", 1.0/minRho);
const colocatedLowerFaceScalarField minRhoInvf("minRhoInv",1.0/ex::interp(minRho));
const colocatedLowerFaceScalarField minRhof("minRhof",ex::interp(minRho));