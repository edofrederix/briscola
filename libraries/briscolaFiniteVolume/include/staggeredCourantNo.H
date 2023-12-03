scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const staggeredScalarLevel& cvl =
        fvMsh.metrics<staggered>().cellVolumes().mshPart();

    staggeredScalarLevel sumPhi(fvMsh,0,Zero);

    forAllCells(sumPhi, d, i, j, k)
        sumPhi(d,i,j,k) = neighborSum(mag(phi(d,i,j,k)));

    CoNum = max(0.5*gMax(sumPhi/cvl)*runTime.deltaTValue());

    meanCoNum =
        max(0.5*(gSum(sumPhi)/gSum(cvl))*runTime.deltaTValue());
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;