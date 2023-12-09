scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const colocatedScalarDirection& cvd =
        fvMsh.metrics<colocated>().cellVolumes().direction();

    colocatedScalarDirection sumPhi(fvMsh,0,0,Zero);

    forAllFaces(sumPhi, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        sumPhi(ijk) += mag(phi(ijk)[fd]);
        sumPhi(nei) += mag(phi(ijk)[fd]);
    }

    CoNum = 0.5*gMax(sumPhi/cvd)*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(cvd))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;