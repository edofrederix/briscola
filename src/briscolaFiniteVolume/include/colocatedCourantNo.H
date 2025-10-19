scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const colocatedScalarDirection& cvd =
        fvMsh.metrics<colocated>().cellVolumes().direction();

    colocatedScalarDirection sumPhi(fvMsh,0,0,Zero);

    forAllFaces(phi, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        sumPhi(ijk) += mag(phi[fd](ijk));
        sumPhi(nei) += mag(phi[fd](ijk));
    }

    CoNum = 0.5*gMax(sumPhi/cvd)*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(cvd))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;