{
    colocatedScalarDirection contError(fvMsh,0,0,Zero);

    const meshField<scalar,colocated>& cv =
        fvMsh.metrics<colocated>().cellVolumes();

    forAllFaces(contError, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        contError(ijk) += phi(ijk)[fd];
        contError(nei) -= phi(ijk)[fd];
    }

    scalar sumLocalContError = 0.0;
    scalar globalContError = 0.0;
    scalar sumVolume = 0.0;

    forAllCells(contError, i, j, k)
    {
        sumLocalContError += Foam::mag(contError(i,j,k))*cv(i,j,k);
        globalContError += contError(i,j,k)*cv(i,j,k);

        sumVolume += cv(i,j,k);
    }

    reduce(sumLocalContError, sumOp<scalar>());
    reduce(globalContError, sumOp<scalar>());
    reduce(sumVolume, sumOp<scalar>());

    sumLocalContError *= runTime.deltaTValue()/sumVolume;
    globalContError *= runTime.deltaTValue()/sumVolume;

    cumulativeContError += globalContError;

    Info<< "time step continuity errors : sum local = " << sumLocalContError
        << ", global = " << globalContError
        << ", cumulative = " << cumulativeContError
        << endl;
}