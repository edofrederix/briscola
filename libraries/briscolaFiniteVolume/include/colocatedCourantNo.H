scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const colocatedScalarDirection& cvd =
        fvMsh.metrics<colocated>().cellVolumes().direction();

    colocatedScalarDirection sumPhi(fvMsh,0,0,Zero);

    forAllCells(sumPhi, i, j, k)
        for (int fd = 0; fd < 3; fd++)
            sumPhi(i,j,k) +=
                mag(phi(i,j,k)[fd])
              + mag(phi(upperNei(i,j,k,fd))[fd]);

    CoNum = 0.5*gMax(sumPhi/cvd)*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(cvd))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;