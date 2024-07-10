const colocatedFaceScalarField& fa =
    fvMsh.metrics<colocated>().faceAreas();

const staggeredLowerFaceScalarField& mu = twoPhase.mu<staggered>();

const staggeredScalarField& rho = twoPhase.rho<staggered>();
