scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const colocatedScalarDirection& cvd =
        fvMsh.metrics<colocated>().cellVolumes().direction();

    colocatedScalarDirection sumPhi(fvMsh,0,0,Zero);

    forAllCells(sumPhi, i, j, k)
        for (int f = 0; f < 6; f++)
            sumPhi(i,j,k) += mag(phi(i,j,k)[f]);

    CoNum = 0.5*gMax(sumPhi/cvd)*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(cvd))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;