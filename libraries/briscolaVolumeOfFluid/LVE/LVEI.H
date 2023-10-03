namespace Foam
{

namespace briscola
{

namespace fv
{

inline scalar LVE::exactCubicSolver
(
    const scalar& Vf,
    const scalar& alpha0,
    const scalar& alpha1,
    const scalar& alpha2,
    const scalar& alpha3,
    const scalar& Cmin,
    const scalar& Cmax,
    const scalar& Vkmin,
    const scalar& Vkmax,
    const scalar& Vtotal
) const
{
    // Analytic cubic equation solver
    //
    // If the leading order term is smaller than a prescribed tolerance or the
    // found root is not in the desired bracket a newton solver is called

    const scalar Cr = Foam::mag(Cmax - Cmin);

    if (Foam::mag(alpha3) > 0)
    {
        if
        (
            Foam::min
            (
                Foam::min
                (
                    Foam::mag(alpha3/alpha0),
                    Foam::mag(alpha3/alpha1)
                ),
                Foam::mag(alpha3/alpha2)
            ) < tol_
        )
        {
            return LVE::newtonCubicSolver
            (
                Vf,
                alpha0,
                alpha1,
                alpha2,
                alpha3,
                Cmin,
                Cmax,
                Vkmin,
                Vkmax,
                Vtotal
            );
        }
        else
        {
            const scalar a0 = (alpha0 - 6.0 * Vf) / alpha3;
            const scalar a1 = alpha1 / alpha3;
            const scalar a2 = alpha2 / alpha3;

            const scalar q = (a1/3) - (Foam::sqr(a2)/9);
            const scalar r = ((a1*a2 - 3*a0)/6) - (Foam::pow3(a2)/27);

            if ((Foam::sqr(r) + Foam::pow3(q)) > 0)
            {
                const scalar A = Foam::cbrt(Foam::mag(r)
                    + Foam::sqrt((Foam::sqr(r) + Foam::pow3(q))));
                scalar t1 = 0;

                if (r >= 0)
                {
                    t1 = A - q/A;
                }
                else
                {
                    t1 = q/A - A;
                }

                const scalar z1 = t1 - a2/3;

                if ((Cmin <= z1 + (tol_ * Cr)) && (z1 <= Cmax + (tol_ * Cr)))
                {
                    if (z1 < Cmin)
                        return -Cmin;
                    else if (z1 > Cmax)
                        return -Cmax;
                    else
                        return -z1;
                }
                else
                {
                    return LVE::newtonCubicSolver
                    (
                        Vf,
                        alpha0,
                        alpha1,
                        alpha2,
                        alpha3,
                        Cmin,
                        Cmax,
                        Vkmin,
                        Vkmax,
                        Vtotal
                    );
                }

            }
            else
            {
                scalar theta = 0;

                if (q < 0)
                {
                    theta = Foam::acos(r/Foam::pow(-q,3/2.));
                }

                const scalar phi1 = theta/3;
                const scalar phi2 =
                    phi1 - ((2*Foam::constant::mathematical::pi)/3);
                const scalar phi3 =
                    phi1 + ((2*Foam::constant::mathematical::pi)/3);

                const scalar z1 = 2 * Foam::sqrt(-q) * Foam::cos(phi1) - a2 /3;
                const scalar z2 = 2 * Foam::sqrt(-q) * Foam::cos(phi2) - a2 /3;
                const scalar z3 = 2 * Foam::sqrt(-q) * Foam::cos(phi3) - a2 /3;

                if ((Cmin <= z1 + (tol_ * Cr)) && (z1 <= Cmax + (tol_ * Cr)))
                {
                    if (z1 < Cmin)
                        return -Cmin;
                    else if (z1 > Cmax)
                        return -Cmax;
                    else
                        return -z1;
                }
                else if
                (
                    (Cmin <= z2 + (tol_ * Cr)) && (z2 <= Cmax + (tol_ * Cr))
                )
                {
                    if (z2 < Cmin)
                        return -Cmin;
                    else if (z2 > Cmax)
                        return -Cmax;
                    else
                        return -z2;
                }
                else if
                (
                    (Cmin <= z3 + (tol_ * Cr)) && (z3 <= Cmax + (tol_ * Cr))
                )
                {
                    if (z3 < Cmin)
                        return -Cmin;
                    else if (z3 > Cmax)
                        return -Cmax;
                    else
                        return -z3;
                }
                else
                {
                    return LVE::newtonCubicSolver
                    (
                        Vf,
                        alpha0,
                        alpha1,
                        alpha2,
                        alpha3,
                        Cmin,
                        Cmax,
                        Vkmin,
                        Vkmax,
                        Vtotal
                    );
                }

            }
        }
    }
    else if (Foam::mag(alpha2) > 0)
    {
        if (Foam::min(Foam::mag(alpha2/alpha0),Foam::mag(alpha2/alpha1)) < tol_)
        {
            return LVE::newtonCubicSolver
            (
                Vf,
                alpha0,
                alpha1,
                alpha2,
                alpha3,
                Cmin,
                Cmax,
                Vkmin,
                Vkmax,
                Vtotal
            );
        }
        else
        {
            const scalar disc = Foam::sqr(alpha1)
                - 4.0 * alpha2 * (alpha0 - 6.0 * Vf);
            if (disc > 0)
            {
                const scalar z1 = (- alpha1 + Foam::sqrt(disc))/(2.0 * alpha2);
                const scalar z2 = (- alpha1 - Foam::sqrt(disc))/(2.0 * alpha2);

                if ((Cmin <= z1 + (tol_ * Cr)) && (z1 <= Cmax + (tol_ * Cr)))
                {
                    if (z1 < Cmin)
                        return -Cmin;
                    else if (z1 > Cmax)
                        return -Cmax;
                    else
                        return -z1;
                }
                else if
                (
                    (Cmin <= z2 + (tol_ * Cr)) && (z2 <= Cmax + (tol_ * Cr))
                )
                {
                    if (z2 < Cmin)
                        return -Cmin;
                    else if (z2 > Cmax)
                        return -Cmax;
                    else
                        return -z2;
                }
                else
                {
                    return LVE::newtonCubicSolver
                    (
                        Vf,
                        alpha0,
                        alpha1,
                        alpha2,
                        alpha3,
                        Cmin,
                        Cmax,
                        Vkmin,
                        Vkmax,
                        Vtotal
                    );
                }
            }
            else
            {
                const scalar z1 = - alpha1 / (2.0 * alpha2);

                if ((Cmin <= z1 + (tol_ * Cr)) && (z1 <= Cmax + (tol_ * Cr)))
                {
                    if (z1 < Cmin)
                        return -Cmin;
                    else if (z1 > Cmax)
                        return -Cmax;
                    else
                        return -z1;
                }
                else
                {
                    return LVE::newtonCubicSolver
                    (
                        Vf,
                        alpha0,
                        alpha1,
                        alpha2,
                        alpha3,
                        Cmin,
                        Cmax,
                        Vkmin,
                        Vkmax,
                        Vtotal
                    );
                }
            }
        }
    }
    else
    {
        const scalar z1 = (6.0 * Vf - alpha0) / alpha1;

        if ((Cmin <= z1 + (tol_ * Cr)) && (z1 <= Cmax + (tol_ * Cr)))
        {
            if (z1 < Cmin)
                return -Cmin;
            else if (z1 > Cmax)
                return -Cmax;
            else
                return -z1;
        }
        else
        {
            return LVE::newtonCubicSolver
            (
                Vf,
                alpha0,
                alpha1,
                alpha2,
                alpha3,
                Cmin,
                Cmax,
                Vkmin,
                Vkmax,
                Vtotal
            );
        }
    }
}

inline scalar LVE::newtonCubicSolver
(
    const scalar& Vf,
    const scalar& alpha0,
    const scalar& alpha1,
    const scalar& alpha2,
    const scalar& alpha3,
    const scalar& Cmin,
    const scalar& Cmax,
    const scalar& Vkmin,
    const scalar& Vkmax,
    const scalar& Vtotal
) const
{
    // Newton method cubic equation solver
    //
    // If the method doesn't converge a bisection solver is called

    scalar Cn = Cmin + (Cmax - Cmin) * (Vf - Vkmax)/(Vkmin - Vkmax);
    scalar fn, dfn;

    for (int iter = 0; iter < 2 * maxIter_; iter++)
    {
        fn = (1.0/6.0) * ((alpha3 * Foam::pow3(Cn))
            + (alpha2 * Foam::sqr(Cn)) + (alpha1 * Cn) + alpha0) - Vf;

        if (Foam::mag(fn/Vtotal) < tol_)
        {
            if (Cn < Cmin)
                return -Cmin;
            else if (Cn > Cmax)
                return -Cmax;
            else
                return -Cn;
            return - Cn;
        }

        dfn = (1.0/6.0) * ((3.0 * alpha3 * Foam::sqr(Cn))
            + (2.0 * alpha2 * Cn) + (alpha1));

        Cn = Cn - fn / (dfn);
    }

    return LVE::bisectionCubicSolver
    (
        Vf,
        alpha0,
        alpha1,
        alpha2,
        alpha3,
        Cmin,
        Cmax,
        Vkmin,
        Vkmax,
        Vtotal
    );
}

inline scalar LVE::bisectionCubicSolver
(
    const scalar& Vf,
    const scalar& alpha0,
    const scalar& alpha1,
    const scalar& alpha2,
    const scalar& alpha3,
    const scalar& Cmin,
    const scalar& Cmax,
    const scalar& Vkmin,
    const scalar& Vkmax,
    const scalar& Vtotal
) const
{
    // Bisection method cubic equation solver

    scalar C1 = Cmin;
    scalar C2 = Cmax;
    scalar Cn = Cmin + (Cmax - Cmin) * (Vf - Vkmax)/(Vkmin - Vkmax);
    scalar fn;

    for (int iter = 0; iter < maxIter_; iter++)
    {
        fn = (1/6.0) * ((alpha3 * Foam::pow3(Cn))
            + (alpha2 * Foam::sqr(Cn)) + (alpha1 * Cn) + alpha0) - Vf;

        if (Foam::mag(fn/Vtotal) < tol_)
        {
            return - Cn;
        }

        if (fn < 0)
        {
            C2 = Cn;
        }
        else
        {
            C1 = Cn;
        }

        Cn = 0.5 * (C1 + C2);
    }

    return - Cn;
}

scalar LVE::cLVE
(
    const vertexVector& v,
    const vector& n,
    const scalar& f,
    const scalar& Vtotal
) const
{
    // CIBRAVE algorithm of J. López, J. Hernández, P. Gómez and F. Faura (2016)
    // adapted for non-convex polyhedrons following, López, J. Hernández, P.
    // Gómez and F. Faura (2019). In this version the cell volume is known.

    scalar Vmin = 0, Vmax;
    scalar Vkmin, Vkmax;
    int Kmin = 0, Kmax = 7, Kf = 0;
    int maxNumberOfFaces = 16;

    scalarList Kjs(maxNumberOfFaces);
    scalarList Mjs(maxNumberOfFaces);
    scalarList Ljs(maxNumberOfFaces);
    int hjs[maxNumberOfFaces];
    int epsilonjs[maxNumberOfFaces];

    vectorList ns(maxNumberOfFaces);
    scalarList Cjs(12);

    int Kmap[8];
    scalarList Cs(8);

    int IA[8];

    // Sorting vertex by their C value

    double originalCs[8];

    for (int i = 0; i < 8; i++)
    {
        originalCs[i] = (n & v[i]);
    }

    double SortedCs[8] = {0};
    int sortedCsIndex[8] = {0};

    int len = 1;
    int place = 0;

    SortedCs[0] = originalCs[0];
    sortedCsIndex[0] = 0;

    for (int i = 1; i < 8; i++)
    {
        place = len;
        for (int j = 0; j < len; j++)
        {
            if (originalCs[i] < SortedCs[j])
            {
                place = j;
                break;
            }
        }

        for (int j = len; j > place; j--)
        {
           SortedCs[j] = SortedCs[j-1];
           sortedCsIndex[j] = sortedCsIndex[j-1];
        }

        SortedCs[place] = originalCs[i];
        sortedCsIndex[place] = i;
        len++;

    }

    // Computation of the normals and their maximum component

    vector aux;

    for (int j = 0; j < 12; j++)
    {
        aux =
            (
                v[faceTriangleDecompCC[j][1]]
              - v[faceTriangleDecompCC[j][0]]
            )
          ^ (
                v[faceTriangleDecompCC[j][2]]
              - v[faceTriangleDecompCC[j][0]]
            );

        ns[j] = aux / Foam::mag(aux);

        Cjs[j] =  - ns[j] & v[faceTriangleDecompCC[j][1]];
    }

    int maxNormalIndex[maxNumberOfFaces];
    scalar auxMax;

    for (int j = 0; j < 12; j++)
    {
        auxMax = Foam::mag(ns[j][0]);
        maxNormalIndex[j] = 0;

        for (int i = 1; i < 3; i++)
        {
            if (Foam::mag(ns[j][i]) > auxMax)
            {
                auxMax = Foam::mag(ns[j][i]);
                maxNormalIndex[j] = i;
            }
        }
    }

    // Check limiting volumes

    Vmax = Vtotal;
    const scalar Vf = Vmax * f;

    if ((Foam::mag(Vf-Vmax)/Vtotal) < tol_)
    {
        return -SortedCs[0];
    }
    else if ((Foam::mag(Vf)/Vtotal) < tol_)
    {
        return -SortedCs[7];
    }

    const scalar Cr = Foam::mag(SortedCs[7] - SortedCs[0]);

    // Prevent edges parallel to the plane from being cut

    Cs[0] = SortedCs[0];
    Kmap[sortedCsIndex[0]] = 0;
    int Knum = 0;

    for (int i = 1; i < 8; i++)
    {
        if (Foam::mag(SortedCs[i-1] - SortedCs[i]) / Cr < tol_)
        {
            Kmap[sortedCsIndex[i]] = Knum;
        }
        else
        {
            Knum++;
            Cs[Knum] = SortedCs[i];
            Kmap[sortedCsIndex[i]] = Knum;
        }
    }

    Kmax = Knum;

    int iterationCount = 0;
    while ((iterationCount < 8) && (Kmin < Kmax))
    {

        iterationCount++;
        int NIPV1[maxNumberOfFaces] = {0};
        int IPV1[maxNumberOfFaces][14];

        int insertedVertex[12][6];
        int totalInsertedVertex = 0;
        int Nf = 12;

        // Find Bracket by interpolation

        int Ip = 8;
        scalar interpolatedC =
            Cs[Kmin] + (Cs[Kmax] - Cs[Kmin]) * (Vf - Vmax)/(Vmin - Vmax);

        for (int i = Kmin; i <= Kmax; i++)
        {
            if (Cs[i] > interpolatedC)
            {
                Kf = i - 1;
                break;
            }
        }

        if (Kf < Kmin)
            Kf = Kmin;
        else if (Kf >= Kmax)
            Kf = Kmax - 1;

        for (int i = 0; i < 8; i++)
        {
            if (Kmap[i] <= Kf)
                IA[i] = 0;
            else
                IA[i] = 1;
        }

        // Arrange new face using Lopez (2020) algorithm

        for (int j = 0; j < 12; j++)
        {
            for (int i = 0; i < 3; i++)
            {
                if (IA[faceTriangleDecompCC[j][i]] == 1)
                    NIPV1[j]++;
            }

            if (NIPV1[j] > 0)
            {
                int index = 0;
                for (int i = 1; i <= 3; i++)
                {
                    int ip1 = faceTriangleDecompCC[j][i % 3];
                    int ip2 = faceTriangleDecompCC[j][(i+1) % 3];

                    if (IA[ip1] == 1)
                    {
                        index++;
                        IPV1[j][index] = ip1;
                    }
                    if (IA[ip1] != IA[ip2])
                    {
                        index++;
                        bool flag = true;

                        for (int aux = 0; aux < totalInsertedVertex; aux ++)
                        {
                            if
                            (
                                ((insertedVertex[aux][1] == ip1)
                                    || (insertedVertex[aux][2] == ip1))
                                & ((insertedVertex[aux][1] == ip2)
                                    || (insertedVertex[aux][2] == ip2))
                            )
                            {
                                flag = false;
                                IPV1[j][index] = insertedVertex[aux][0];
                                if (IA[ip2] == 1)
                                {
                                    insertedVertex[aux][3] = j;
                                    insertedVertex[aux][4] = index;
                                }
                                break;
                            }
                        }

                        if (flag)
                        {
                            insertedVertex[totalInsertedVertex][0] = Ip;
                            IPV1[j][index] = Ip;

                            if (IA[ip2] == 1)
                            {
                                insertedVertex[totalInsertedVertex][3] = j;
                                insertedVertex[totalInsertedVertex][4] = index;
                                insertedVertex[totalInsertedVertex][2] = ip1;
                                insertedVertex[totalInsertedVertex][1] = ip2;
                            }
                            else
                            {
                                insertedVertex[totalInsertedVertex][1] = ip1;
                                insertedVertex[totalInsertedVertex][2] = ip2;
                            }

                            insertedVertex[totalInsertedVertex][5] = 0;
                            totalInsertedVertex++;
                            Ip++;
                        }
                    }
                }

                NIPV1[j] = index;
                IPV1[j][0] = IPV1[j][NIPV1[j]];
                IPV1[j][NIPV1[j]+1] = IPV1[j][1];
            }
        }

        int index = 1;
        bool untaggedVertex = true;
        int taggedVertex = 1;

        IPV1[Nf][index] = 8;
        insertedVertex[0][5] = 1;
        int counter = 0;

        while (untaggedVertex)
        {
            if (counter > 13)
            {
                FatalErrorInFunction
                    << " Arrangement of the new face in LVE failed."
                    << endl << abort(FatalError);
            }

            int auxj = insertedVertex[IPV1[Nf][index]-8][3];
            int auxi = insertedVertex[IPV1[Nf][index]-8][4];

            if (IPV1[auxj][auxi - 1] != IPV1[Nf][1])
            {
                index++;
                IPV1[Nf][index] = IPV1[auxj][auxi - 1];
                insertedVertex[IPV1[Nf][index]-8][5] = 1;
                taggedVertex++;
            }
            else
            {
                NIPV1[Nf] = index;
                IPV1[Nf][0] = IPV1[Nf][NIPV1[Nf]];
                IPV1[Nf][NIPV1[Nf]+1] = IPV1[Nf][1];
                Nf++;

                if (taggedVertex == totalInsertedVertex)
                {
                    untaggedVertex = false;
                }
                else
                {
                    index = 1;
                    for (int auxip = 0; auxip < totalInsertedVertex; auxip++)
                    {
                        if (insertedVertex[auxip][5] == 0)
                        {
                            IPV1[Nf][index] = insertedVertex[auxip][0];
                            insertedVertex[auxip][5] = 1;
                            taggedVertex++;
                            break;
                        }
                    }
                }

            }

            counter++;
        }

        // Compute the piecewise analytical function

        scalar alpha0 = 0;
        scalar alpha1 = 0;
        scalar alpha2 = 0;
        scalar alpha3 = 0;

        for (int j = 12; j < Nf; j++)
        {
            ns[j] = - n;
        }

        vectorList x0(Ip);
        vectorList E(Ip);

        for (int i = 0; i < Ip; i++)
        {
            if (i < 8)
            {
                x0[i] = v[i];
                E[i] = Foam::zero();
            }
            else
            {
                aux = v[insertedVertex[i-8][2]]
                    - v[insertedVertex[i-8][1]];
                aux = aux / Foam::mag(aux);
                E[i] = aux / (n & aux);
                x0[i] = v[insertedVertex[i-8][1]]
                    - (n & v[insertedVertex[i-8][1]]) * E[i];
            }
        }

        for (int j = 12; j < Nf; j++)
        {
            auxMax = Foam::mag(ns[j][0]);
            maxNormalIndex[j] = 0;

            for (int i = 1; i < 3; i++)
            {
                if (Foam::mag(ns[j][i]) > auxMax)
                {
                    auxMax = Foam::mag(ns[j][i]);
                    maxNormalIndex[j] = i;
                }
            }
        }

        for (int j = 0; j < Nf; j++)
        {
            if (NIPV1[j] > 0)
            {
                hjs[j] = int((NIPV1[j]-2)/2);
                epsilonjs[j] = NIPV1[j] % 2;

                Kjs[j] = Foam::zero();
                Ljs[j] = Foam::zero();
                Mjs[j] = Foam::zero();

                for (int i = 2; i <= (hjs[j] + 1); i++)
                {
                    Kjs[j] +=
                        crossProductComponent
                        (
                            (x0[IPV1[j][2*i - 1]] - x0[IPV1[j][1]]),
                            (x0[IPV1[j][2*i]] - x0[IPV1[j][2*i-2]]),
                            maxNormalIndex[j]
                        );

                    Ljs[j] +=
                        crossProductComponent
                        (
                            (x0[IPV1[j][2*i - 1]] - x0[IPV1[j][1]]),
                            (E[IPV1[j][2*i]] - E[IPV1[j][2*i-2]]),
                            maxNormalIndex[j]
                        )
                      + crossProductComponent
                        (
                            (E[IPV1[j][2*i - 1]] - E[IPV1[j][1]]),
                            (x0[IPV1[j][2*i]] - x0[IPV1[j][2*i-2]]),
                            maxNormalIndex[j]
                        );

                    Mjs[j] +=
                        crossProductComponent
                        (
                            (E[IPV1[j][2*i - 1]] - E[IPV1[j][1]]),
                            (E[IPV1[j][2*i]] - E[IPV1[j][2*i-2]]),
                            maxNormalIndex[j]
                        );
                }

                if (epsilonjs[j] == 1)
                {
                    Kjs[j] +=
                        crossProductComponent
                        (
                            (x0[IPV1[j][NIPV1[j]]] - x0[IPV1[j][1]]),
                            (x0[IPV1[j][1]] - x0[IPV1[j][NIPV1[j]-1]]),
                            maxNormalIndex[j]
                        );

                    Ljs[j] +=
                        crossProductComponent
                        (
                            (x0[IPV1[j][NIPV1[j]]] - x0[IPV1[j][1]]),
                            (E[IPV1[j][1]] - E[IPV1[j][NIPV1[j]-1]]),
                            maxNormalIndex[j]
                        )
                      + crossProductComponent
                        (
                            (E[IPV1[j][NIPV1[j]]] - E[IPV1[j][1]]),
                            (x0[IPV1[j][1]] - x0[IPV1[j][NIPV1[j]-1]]),
                            maxNormalIndex[j]
                        );

                    Mjs[j] +=
                        crossProductComponent
                        (
                            (E[IPV1[j][NIPV1[j]]] - E[IPV1[j][1]]),
                            (E[IPV1[j][1]] - E[IPV1[j][NIPV1[j]-1]]),
                            maxNormalIndex[j]
                        );
                }
            }
        }

        for (int j = 0; j < 12; j++)
        {
            if (NIPV1[j] > 0)
            {
                alpha0 += - Cjs[j] * (Kjs[j] / ns[j][maxNormalIndex[j]]);
                alpha1 += - Cjs[j] * (Ljs[j] / ns[j][maxNormalIndex[j]]);
                alpha2 += - Cjs[j] * (Mjs[j] / ns[j][maxNormalIndex[j]]);
            }
        }

        for (int j = 12; j < Nf; j++)
        {
            if (NIPV1[j] > 0)
            {
                alpha1 += - (Kjs[j] / ns[j][maxNormalIndex[j]]);
                alpha2 += - (Ljs[j] / ns[j][maxNormalIndex[j]]);
                alpha3 += - (Mjs[j] / ns[j][maxNormalIndex[j]]);
            }
        }

        // Check bracket and solve inverse problem

        Vkmin = (1/6.0) * ((alpha3 * Foam::pow3(Cs[Kf + 1]))
            + (alpha2 * Foam::sqr(Cs[Kf + 1])) + (alpha1 * Cs[Kf + 1]) + alpha0);
        Vkmax = (1/6.0) * ((alpha3 * Foam::pow3(Cs[Kf]))
            + (alpha2 * Foam::sqr(Cs[Kf])) + (alpha1 * Cs[Kf]) + alpha0);

        if (Foam::mag(Vkmin - Vf)/Vtotal < tol_)
        {
            return -Cs[Kf + 1];
        }
        else if (Foam::mag(Vkmax - Vf)/Vtotal < tol_)
        {
            return -Cs[Kf];
        }
        else if ((Vkmin <= Vf) & (Vf <= Vkmax))
        {

            // The equation is solved using an exact solver except in some cases
            // where the leading order coefficient is small. In those cases the
            // newton method is used and if the newton method doesn't converge
            // the bisection method is used

            scalar z = LVE::exactCubicSolver
            (
                Vf,
                alpha0,
                alpha1,
                alpha2,
                alpha3,
                Cs[Kf],
                Cs[Kf + 1],
                Vkmin,
                Vkmax,
                Vtotal
            );

            if (z != z)
            {
                FatalErrorInFunction
                    << "LVE solver failed. Cn: " << z
                    << endl << abort(FatalError);
            }

            return z;
        }
        else if (Vf > Vkmax)
        {
            if (Kmin == Kf)
            {
                return -Cs[Kf];
            }
            else
            {
                Vmin = Vkmax;
                Kmax = Kf;
            }
        }
        else
        {
            if (Kmax == (Kf + 1))
            {
                return -Cs[Kf + 1];
            }
            else
            {
                Vmax = Vkmin;
                Kmin = Kf + 1;
            }
        }

    }

    FatalErrorInFunction
        <<  "LVE solver failed." << endl << abort(FatalError);

    return Cs[0];
}

inline scalar LVE::cLVE
(
    const vertexVector& v,
    const vector& n,
    const scalar& f
) const
{

    // Compute the volume of a general cell and then call the CIBRAVE algorithm

    scalarList Ajs(12);
    int hjs[12];
    int epsilonjs[12];

    vectorList ns(12);
    scalarList Cjs(12);
    vector aux;

    for (int j = 0; j < 12; j++)
    {
        aux =
            (
                v[faceTriangleDecompCC[j][1]]
              - v[faceTriangleDecompCC[j][0]]
            )
          ^ (
                v[faceTriangleDecompCC[j][2]]
              - v[faceTriangleDecompCC[j][0]]
            );

        ns[j] = aux / Foam::mag(aux);

        Cjs[j] =  - ns[j] & v[faceTriangleDecompCC[j][0]];
    }

    int maxNormalIndex[12];
    scalar auxMax;

    for (int j = 0; j < 12; j++)
    {
        auxMax = Foam::mag(ns[j][0]);
        maxNormalIndex[j] = 0;

        for (int i = 1; i < 3; i++)
        {
            if (Foam::mag(ns[j][i]) > auxMax)
            {
                auxMax = Foam::mag(ns[j][i]);
                maxNormalIndex[j] = i;
            }
        }
    }

    for (int j = 0; j < 12; j++)
    {
        hjs[j] = int((3-2)/2);
        epsilonjs[j] = 3 % 2;

        Ajs[j] = Foam::zero();

        for (int i = 1; i < (hjs[j] + 1); i++)
        {
            Ajs[j] +=
                crossProductComponent
                (
                    (
                        v[faceTriangleDecompCC[j][(2*i - 1) % 3]]
                      - v[faceTriangleDecompCC[j][0]]
                    ),
                    (
                        v[faceTriangleDecompCC[j][(2*i) % 3]]
                      - v[faceTriangleDecompCC[j][(2*i-2) % 3]]
                    ),
                    maxNormalIndex[j]
                );
        }

        if (epsilonjs[j] == 1)
        {
            Ajs[j] +=
                crossProductComponent
                (
                    (
                        v[faceTriangleDecompCC[j][2]] -
                        v[faceTriangleDecompCC[j][0]]
                    ),
                    (
                        v[faceTriangleDecompCC[j][0]] -
                        v[faceTriangleDecompCC[j][1]]
                    ),
                    maxNormalIndex[j]
                );
        }
    }

    scalar Vtotal = 0;
    for (int j = 0; j < 12; j++)
    {
        Vtotal += - Cjs[j] * (Ajs[j] / ns[j][maxNormalIndex[j]]);
    }

    Vtotal /= 6;

    return cLVE(v,n,f,Vtotal);
}

inline scalar LVE::pLVE
(
    const vertexVector& vertices,
    const vector& no,
    const scalar& fo
) const
{
    // Algorithm of Scardovelli & Zaleski (2000) as implemented in Basilisk and
    // VOFTools

    #ifdef FULLDEBUG

    vector d
    (
        Foam::mag(vertices.rba() - vertices.lba()),
        Foam::mag(vertices.lta() - vertices.lba()),
        Foam::mag(vertices.lbf() - vertices.lba())
    );

    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
    {
        const label v0 = vertexNumsInEdge[i*4+j][0];
        const label v1 = vertexNumsInEdge[i*4+j][1];
        const scalar dj = Foam::mag(vertices[v0] - vertices[v1]);

        if (Foam::mag(dj - d[i]) > 1e-12)
            FatalErrorInFunction
                << "Cell is not a parallelepiped" << endl
                << abort(FatalError);
    }

    #endif

    const tensor T
    (
        vertices.rba() - vertices.lba(),
        vertices.lta() - vertices.lba(),
        vertices.lbf() - vertices.lba()
    );

    vector n = T & no;
    const scalar S = Foam::mag(n);
    n = n/S;

    scalarList Cs(8);

    for (int i = 0; i < 8; i++)
    {
        Cs[i] = -(n & unitCube[i]);
    }

    const scalar CMin = Foam::min(Cs);
    const scalar CMax = Foam::max(Cs);

    // Solve inverse if f > 0.5

    const scalar f = Foam::max(Foam::min(fo, 1.0-fo),0);

    scalar s = cmptSum(cmptMag(n));
    vector m = cmptMag(n)/s;

    scalarList ml(3);

    for (int i = 0; i < 3; i++)
        ml[i] = m[i];

    sort(ml);

    const scalar& m1 = ml[0];
    const scalar& m2 = ml[1];
    const scalar& m3 = ml[2];

    const scalar m12 = m1 + m2;

    const scalar pr = Foam::max(6.0*m1*m2*m3, 1e-50);
    const scalar v1 = Foam::pow3(m1)/pr;
    const scalar v2 = v1 + (m2 - m1)/(2.0*m3);

    scalar v3, mm;

    if (m3 < m12)
    {
        mm = m3;
        v3 =
            (
              + (3*m12 - m3)*Foam::sqr(m3)
              + (m1 - 3*m3 )*Foam::sqr(m1)
              + (m2 - 3*m3 )*Foam::sqr(m2)
            )
          / pr;
    }
    else
    {
        mm = m12;
        v3 = m12/(2.0*m3);
    }

    scalar alpha = 0;

    if (f < v1)
    {
        alpha = Foam::pow(pr*f, 1.0/3.0);
    }
    else if (f < v2)
    {
        alpha = 0.5 * (m1 + Foam::sqrt(Foam::sqr(m1) + 8.0*m2*m3*(f - v1)));
    }
    else if (f < v3)
    {
        const scalar p12 = Foam::sqrt(2.0*m1*m2);
        const scalar q = 3.0*(m12 - 2.0*m3*f)/(4.0*p12);
        const scalar theta = Foam::acos(Foam::max(Foam::min(q,1),-1))/3.0;
        const scalar cs = Foam::cos(theta);

        alpha = p12*(Foam::sqrt(3.0*(1.0 - Foam::sqr(cs))) - cs) + m12;
    }
    else if (m12 <= m3)
    {
        alpha = m3*f + 0.5*mm;
    }
    else
    {
        const scalar p = m1*(m2 + m3) + m2*m3 - 0.25;
        const scalar p12 = Foam::sqrt(p);
        const scalar q = 3.0*m1*m2*m3*(0.5 - f)/(2.0*p*p12);
        const scalar theta = Foam::acos(Foam::max(Foam::min(q,1),-1))/3.0;
        const scalar cs = Foam::cos(theta);

        alpha = p12*(Foam::sqrt(3.0*(1.0 - Foam::sqr(cs))) - cs) + 0.5;
    }

    scalar C;

    if (fo <= 0.5)
    {
        C = CMin + alpha * Foam::mag(CMax - CMin);
    }
    else
    {
        C = CMax - alpha * Foam::mag(CMax - CMin);
    }

    return C*S - (no & vertices.lba());
}

}

}

}
