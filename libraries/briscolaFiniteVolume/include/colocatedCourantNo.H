scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const colocatedScalarDirection& cvd =
        fvMsh.metrics<colocated>().cellVolumes().direction();

    colocatedScalarDirection sumPhi(fvMsh,0,0);

    forAllCells(sumPhi, i, j, k)
        sumPhi(i,j,k) = neighborSum(mag(phi(i,j,k)));

    CoNum = 0.5*gMax(sumPhi/cvd)*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(cvd))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;