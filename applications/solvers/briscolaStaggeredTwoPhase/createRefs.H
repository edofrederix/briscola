const colocatedFaceScalarField& fa =
    fvMsh.metrics<colocated>().faceAreas();

const staggeredLowerFaceScalarField& mu = twoPhase.mu<staggered>();
