{
    scalarBlock contError(fvMsh.N());

    const meshDirection<scalar,colocated>& cv =
        fvMsh.metrics<colocated>().cellVolumes()[0][0];

    const meshDirection<hexScalar,colocated>& phid = phi[0][0];

    forAllBlock(contError, i, j, k)
    {
        contError(i,j,k) = neighborSkewSum(phid(i,j,k));
    }

    scalar sumLocalContError = 0.0;
    scalar globalContError = 0.0;
    scalar sumVolume = 0.0;

    forAllBlock(contError, i, j, k)
    {
        sumLocalContError += Foam::mag(contError(i,j,k))*cv(i,j,k);
        globalContError += contError(i,j,k)*cv(i,j,k);

        sumVolume += cv(i,j,k);
    }

    sumLocalContError *= runTime.deltaTValue()/sumVolume;
    globalContError *= runTime.deltaTValue()/sumVolume;

    cumulativeContError += globalContError;

    Info<< "time step continuity errors : sum local = " << sumLocalContError
        << ", global = " << globalContError
        << ", cumulative = " << cumulativeContError
        << endl;
}