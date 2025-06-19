namespace Foam
{

namespace briscola
{

namespace fv
{

inline scalar LVE::exactCubicSolver
(
    const scalar Vf,
    const scalarList& alpha,
    const scalar Cmin,
    const scalar Cmax,
    const scalar Wmin,
    const scalar Wmax,
    const scalar V
) const
{
    // Analytic cubic equation solver. If the leading order term is smaller than
    // a prescribed tolerance or the found root is not in the desired bracket a
    // Newton solver is called

    const scalar dC = Foam::mag(Cmax - Cmin);

    if (Foam::mag(alpha[3]) > 0)
    {
        if
        (
            Foam::min
            (
                Foam::min
                (
                    Foam::mag(alpha[3]/stabilise(alpha[0],1e-16)),
                    Foam::mag(alpha[3]/stabilise(alpha[1],1e-16))
                ),
                Foam::mag(alpha[3]/stabilise(alpha[2],1e-16))
            ) < tol_
        )
        {
            return LVE::NewtonCubicSolver
            (
                Vf,
                alpha,
                Cmin,
                Cmax,
                Wmin,
                Wmax,
                V
            );
        }
        else
        {
            const scalar a0 = (alpha[0] - 6.0*Vf)/alpha[3];
            const scalar a1 = alpha[1]/alpha[3];
            const scalar a2 = alpha[2]/alpha[3];

            const scalar q = a1/3 - Foam::sqr(a2)/9.0;
            const scalar r = (a1*a2 - 3*a0)/6.0 - Foam::pow3(a2)/27.0;

            if (Foam::sqr(r) + Foam::pow3(q) > 0)
            {
                const scalar A =
                    Foam::cbrt
                    (
                        Foam::mag(r)
                      + Foam::sqrt(Foam::sqr(r) + Foam::pow3(q))
                    );

                const scalar t1 = r >= 0 ? A - q/A : q/A - A;
                const scalar z1 = t1 - a2/3;

                if (Cmin <= z1 + tol_*dC && z1 <= Cmax + tol_*dC)
                {
                    return -Foam::max(Foam::min(z1,Cmax),Cmin);
                }
                else
                {
                    return LVE::NewtonCubicSolver
                    (
                        Vf,
                        alpha,
                        Cmin,
                        Cmax,
                        Wmin,
                        Wmax,
                        V
                    );
                }
            }
            else
            {
                const scalar theta =
                    q < 0 ? Foam::acos(r/Foam::pow(-q,3.0/2.0)) : 0.0;

                const scalar pi = Foam::constant::mathematical::pi;

                const scalar phi1 = theta/3.0;
                const scalar phi2 = phi1 - 2.0*pi/3.0;
                const scalar phi3 = phi1 + 2.0*pi/3.0;

                const scalar z1 = 2.0*Foam::sqrt(-q)*Foam::cos(phi1) - a2/3.0;
                const scalar z2 = 2.0*Foam::sqrt(-q)*Foam::cos(phi2) - a2/3.0;
                const scalar z3 = 2.0*Foam::sqrt(-q)*Foam::cos(phi3) - a2/3.0;

                if (Cmin <= z1 + tol_*dC && z1 <= Cmax + tol_*dC)
                {
                    return -Foam::max(Foam::min(z1,Cmax),Cmin);
                }
                else if (Cmin <= z2 + tol_*dC && z2 <= Cmax + tol_*dC)
                {
                    return -Foam::max(Foam::min(z2,Cmax),Cmin);
                }
                else if (Cmin <= z3 + tol_*dC && z3 <= Cmax + tol_*dC)
                {
                    return -Foam::max(Foam::min(z3,Cmax),Cmin);
                }
                else
                {
                    return LVE::NewtonCubicSolver
                    (
                        Vf,
                        alpha,
                        Cmin,
                        Cmax,
                        Wmin,
                        Wmax,
                        V
                    );
                }
            }
        }
    }
    else if (Foam::mag(alpha[2]) > 0)
    {
        if
        (
            Foam::min
            (
                Foam::mag(alpha[2]/stabilise(alpha[0],1e-16)),
                Foam::mag(alpha[2]/stabilise(alpha[1],1e-16))
            ) < tol_
        )
        {
            return LVE::NewtonCubicSolver
            (
                Vf,
                alpha,
                Cmin,
                Cmax,
                Wmin,
                Wmax,
                V
            );
        }
        else
        {
            const scalar disc =
                Foam::sqr(alpha[1]) - 4.0*alpha[2]*(alpha[0] - 6.0*Vf);

            if (disc > 0)
            {
                const scalar z1 = 0.5*(-alpha[1] + Foam::sqrt(disc))/alpha[2];
                const scalar z2 = 0.5*(-alpha[1] - Foam::sqrt(disc))/alpha[2];

                if (Cmin <= z1 + tol_*dC && z1 <= Cmax + tol_*dC)
                {
                    return -Foam::max(Foam::min(z1,Cmax),Cmin);
                }
                else if (Cmin <= z2 + tol_*dC && z2 <= Cmax + tol_*dC)
                {
                    return -Foam::max(Foam::min(z2,Cmax),Cmin);
                }
                else
                {
                    return LVE::NewtonCubicSolver
                    (
                        Vf,
                        alpha,
                        Cmin,
                        Cmax,
                        Wmin,
                        Wmax,
                        V
                    );
                }
            }
            else
            {
                const scalar z1 = -0.5*alpha[1]/alpha[2];

                if (Cmin <= z1 + tol_*dC && z1 <= Cmax + tol_*dC)
                {
                    return -Foam::max(Foam::min(z1,Cmax),Cmin);
                }
                else
                {
                    return LVE::NewtonCubicSolver
                    (
                        Vf,
                        alpha,
                        Cmin,
                        Cmax,
                        Wmin,
                        Wmax,
                        V
                    );
                }
            }
        }
    }
    else
    {
        const scalar z1 = (6.0*Vf - alpha[0])/alpha[1];

        if (Cmin <= z1 + tol_*dC && z1 <= Cmax + tol_*dC)
        {
            return -Foam::max(Foam::min(z1,Cmax),Cmin);
        }
        else
        {
            return LVE::NewtonCubicSolver
            (
                Vf,
                alpha,
                Cmin,
                Cmax,
                Wmin,
                Wmax,
                V
            );
        }
    }
}

inline scalar LVE::NewtonCubicSolver
(
    const scalar Vf,
    const scalarList& alpha,
    const scalar Cmin,
    const scalar Cmax,
    const scalar Wmin,
    const scalar Wmax,
    const scalar V
) const
{
    // Newton method cubic equation solver. If the method doesn't converge a
    // bisection solver is called.

    scalar C = Cmin + (Cmax - Cmin)*(Wmax - Vf)/(Wmax - Wmin);

    for (label iter = 0; iter < 2*maxIter_; iter++)
    {
        const scalar fn =
            (
                alpha[3]*Foam::pow3(C)
              + alpha[2]*Foam::sqr(C)
              + alpha[1]*C
              + alpha[0]
            )/6.0
          - Vf;

        if (Foam::mag(fn/V) < tol_)
            return -Foam::max(Foam::min(C,Cmax),Cmin);

        const scalar dfn =
            (
                alpha[3]*Foam::sqr(C)*3.0
              + alpha[2]*C*2.0
              + alpha[1]
            )/6.0;

        C = C - fn/stabilise(dfn,1e-25);
    }

    return LVE::bisectionCubicSolver
    (
        Vf,
        alpha,
        Cmin,
        Cmax,
        Wmin,
        Wmax,
        V
    );
}

inline scalar LVE::bisectionCubicSolver
(
    const scalar Vf,
    const scalarList& alpha,
    const scalar Cmin,
    const scalar Cmax,
    const scalar Wmin,
    const scalar Wmax,
    const scalar V
) const
{
    // Bisection method cubic equation solver. Last resort.

    scalar C1 = Cmin;
    scalar C2 = Cmax;
    scalar C = Cmin + (Cmax - Cmin)*(Wmax - Vf)/(Wmax - Wmin);

    for (label iter = 0; iter < maxIter_; iter++)
    {
        const scalar fn =
            (
                alpha[3]*Foam::pow3(C)
              + alpha[2]*Foam::sqr(C)
              + alpha[1]*C
              + alpha[0]
            )/6.0
          - Vf;

        if (Foam::mag(fn/V) < tol_)
            return -Foam::max(Foam::min(C,Cmax),Cmin);

        // Adjust limits

        if (fn < 0)
        {
            C2 = C;
        }
        else
        {
            C1 = C;
        }

        // Bisection

        C = 0.5*(C1 + C2);
    }

    FatalError
        << "Could not find LVE solution" << endl << abort(FatalError);

    return -C;
}

inline scalar LVE::cLVE
(
    const vertexVector& v,
    const vector& n,
    const scalar f,
    const scalar V,
    const vectorList& triAreas,
    const vectorList& triNormals,
    const scalarList& Cj,
    const labelList& m
) const
{
    // CIBRAVE algorithm of J. López, J. Hernández, P. Gómez and F. Faura (2016)
    // adapted for non-convex polyhedrons following, López, J. Hernández, P.
    // Gómez and F. Faura (2019). In this version the cell volume is known.

    const scalar Vf = V*f;

    // Sorting vertex by their C value. The clean function trims truncation
    // error.

    scalarList C(8);
    forAll(C, i)
        C[i] = trimPrecision(n & v[i]);

    labelList order;
    sortedOrder(C, order);

    scalarList Cs(8);
    forAll(C, i)
        Cs[i] = C[order[i]];

    // Check limiting volumes

    if (f < tol_)
    {
        return -Cs.last();
    }
    else if (f > 1.0-tol_)
    {
        return -Cs.first();
    }

    // Unique sorting

    labelList uniqueOrder(8, 0);
    for (label i = 1; i < 8; i++)
        uniqueOrder[i] =
            uniqueOrder[i-1] + (Cs[i] != Cs[i-1] ? 1 : 0);

    // Unique C values

    scalarList Cu(uniqueOrder.last()+1);
    forAll(Cu, i)
        Cu[i] = Cs[findIndex(uniqueOrder,i)];

    label Kmin = 0, Kmax = Cu.size()-1;
    scalar Vmin = 0.0, Vmax = V;

    label iter = 0;
    while (iter < 8 && Kmin < Kmax)
    {
        iter++;

        // Find bracket by interpolation

        scalar Ci =
            Cu[Kmin] + (Cu[Kmax] - Cu[Kmin])*(Vmax - Vf)/(Vmax - Vmin);

        // Lower unique vertex index

        const label Kj =
            Foam::min(Foam::max(findLower(Cu,Ci), Kmin), Kmax-1);

        // Check once every iteration if vertex C values are larger than Ci

        List<bool> IA(8);
        forAll(IA, i)
            IA[i] = C[i] > Ci;

        // Arrange new face using Lopez (2020) algorithm. The maximum number of
        // faces is 16.

        labelList NIPV1(16, 0);
        List<labelList> IPV1(16, labelList(14,-1));
        List<labelList> IN(0);

        // Total number of vertices

        label I = 8;

        // Loop over face triangles

        for (int j = 0; j < 12; j++)
        {
            // In this face triangle, count the number of vertices that have a C
            // value larger than Ci

            for (int i = 0; i < 3; i++)
                NIPV1[j] += IA[faceTriangleDecompCC[j][i]];

            // Check if at least one vertex has a large C value

            if (NIPV1[j] > 0)
            {
                label k = 0;
                for (label i = 0; i < 3; i++)
                {
                    label i1 = faceTriangleDecompCC[j][(i+1) % 3];
                    label i2 = faceTriangleDecompCC[j][(i+2) % 3];

                    if (IA[i1])
                    {
                        k++;
                        IPV1[j][k] = i1;
                    }

                    if (IA[i1] != IA[i2])
                    {
                        k++;
                        bool add = true;

                        forAll(IN, l)
                        {
                            if
                            (
                                (IN[l][1] == i1 || IN[l][2] == i1)
                             && (IN[l][1] == i2 || IN[l][2] == i2)
                            )
                            {
                                add = false;

                                IPV1[j][k] = IN[l][0];

                                if (IA[i2] == 1)
                                {
                                    IN[l][3] = j;
                                    IN[l][4] = k;
                                }

                                break;
                            }
                        }

                        if (add)
                        {
                            IN.append(labelList(6));
                            IN.last()[0] = I;
                            IPV1[j][k] = I;

                            if (IA[i2] == 1)
                            {
                                IN.last()[1] = i2;
                                IN.last()[2] = i1;
                                IN.last()[3] = j;
                                IN.last()[4] = k;
                            }
                            else
                            {
                                IN.last()[1] = i1;
                                IN.last()[2] = i2;
                                IN.last()[3] = -1;
                                IN.last()[4] = -1;
                            }

                            IN.last()[5] = 0;
                            I++;
                        }
                    }
                }

                NIPV1[j] = k;
                IPV1[j][0] = IPV1[j][NIPV1[j]];
                IPV1[j][NIPV1[j]+1] = IPV1[j][1];
            }
        }

        // Total number of faces

        label J = 12;

        // Arrange faces

        label K = 1, nTags = 1;
        IPV1[J][K] = 8;
        IN[0][5] = 1;

        label count = 0;
        bool keepGoing = true;
        while (keepGoing)
        {
            if (count > 13)
                FatalErrorInFunction
                    << "Arrangement of the new face in LVE failed."
                    << endl << abort(FatalError);

            int i = IN[IPV1[J][K]-8][4];
            int j = IN[IPV1[J][K]-8][3];

            if (IPV1[j][i-1] != IPV1[J][1])
            {
                K++;

                IPV1[J][K] = IPV1[j][i-1];
                IN[IPV1[J][K]-8][5] = 1;
                nTags++;
            }
            else
            {
                NIPV1[J] = K;
                IPV1[J][0] = IPV1[J][NIPV1[J]];
                IPV1[J][NIPV1[J]+1] = IPV1[J][1];

                J++;

                if (nTags == IN.size())
                {
                    keepGoing = false;
                }
                else
                {
                    K = 1;
                    for (int k = 0; k < IN.size(); k++)
                    {
                        if (IN[k][5] == 0)
                        {
                            IPV1[J][K] = IN[k][0];
                            IN[k][5] = 1;

                            nTags++;

                            break;
                        }
                    }
                }
            }

            count++;
        }

        // Extend list of normals

        vectorList normals(triNormals);
        normals.resize(16);

        // Set new face normals to the interface normal

        for (int j = 12; j < J; j++)
            normals[j] = -n;

        // Extend list of vertices

        vectorList w(I);
        vectorList E(I,Zero);

        for (label i = 0; i < 8; i++)
            w[i] = v[i];

        for (int i = 8; i < I; i++)
        {
            const vector z = v[IN[i-8][2]] - v[IN[i-8][1]];
            const vector normal = z/Foam::mag(z);

            E[i] = normal/(n & normal);
            w[i] = v[IN[i-8][1]] - (n & v[IN[i-8][1]]) * E[i];
        }

        // Extend list of m indices

        labelList M(m);
        M.resize(J);

        for (int j = 12; j < J; j++)
            M[j] =
                Foam::mag(normals[j].y()) > Foam::mag(normals[j].x())
              ? Foam::mag(normals[j].z()) > Foam::mag(normals[j].y()) ? 2 : 1
              : Foam::mag(normals[j].z()) > Foam::mag(normals[j].x()) ? 2 : 0;

        // Set area vector components

        List<scalarList> areas(3, scalarList(16, 0.0));

        for (int j = 0; j < J; j++)
        {
            if (NIPV1[j] > 0)
            {
                for (int i = 2; i <= (label((NIPV1[j]-2)/2) + 1); i++)
                {
                    areas[0][j] +=
                        crossProductComponent
                        (
                            (w[IPV1[j][2*i-1]] - w[IPV1[j][1    ]]),
                            (w[IPV1[j][2*i  ]] - w[IPV1[j][2*i-2]]),
                            M[j]
                        );

                    areas[1][j] +=
                        crossProductComponent
                        (
                            (w[IPV1[j][2*i - 1]] - w[IPV1[j][1    ]]),
                            (E[IPV1[j][2*i    ]] - E[IPV1[j][2*i-2]]),
                            M[j]
                        )
                      + crossProductComponent
                        (
                            (E[IPV1[j][2*i - 1]] - E[IPV1[j][1    ]]),
                            (w[IPV1[j][2*i    ]] - w[IPV1[j][2*i-2]]),
                            M[j]
                        );

                    areas[2][j] +=
                        crossProductComponent
                        (
                            (E[IPV1[j][2*i - 1]] - E[IPV1[j][1    ]]),
                            (E[IPV1[j][2*i    ]] - E[IPV1[j][2*i-2]]),
                            M[j]
                        );
                }

                if (NIPV1[j] % 2)
                {
                    areas[0][j] +=
                        crossProductComponent
                        (
                            (w[IPV1[j][NIPV1[j]]] - w[IPV1[j][1         ]]),
                            (w[IPV1[j][1       ]] - w[IPV1[j][NIPV1[j]-1]]),
                            M[j]
                        );

                    areas[1][j] +=
                        crossProductComponent
                        (
                            (w[IPV1[j][NIPV1[j]]] - w[IPV1[j][1         ]]),
                            (E[IPV1[j][1       ]] - E[IPV1[j][NIPV1[j]-1]]),
                            M[j]
                        )
                      + crossProductComponent
                        (
                            (E[IPV1[j][NIPV1[j]]] - E[IPV1[j][1         ]]),
                            (w[IPV1[j][1       ]] - w[IPV1[j][NIPV1[j]-1]]),
                            M[j]
                        );

                    areas[2][j] +=
                        crossProductComponent
                        (
                            (E[IPV1[j][NIPV1[j]]] - E[IPV1[j][1         ]]),
                            (E[IPV1[j][1       ]] - E[IPV1[j][NIPV1[j]-1]]),
                            M[j]
                        );
                }
            }
        }

        // Alpha values

        scalarList alpha(4, 0.0);

        for (int j = 0; j < 12; j++)
            if (NIPV1[j] > 0)
                for (int i = 0; i < 3; i++)
                    alpha[i] -=
                        Cj[j]*areas[i][j]/normals[j][M[j]];

        for (int j = 12; j < J; j++)
            if (NIPV1[j] > 0)
                for (int i = 0; i < 3; i++)
                    alpha[i+1] -=
                        areas[i][j]/normals[j][M[j]];

        // Check bracket volumes and solve inverse problem

        const scalar Wmin =
            (
                alpha[3]*Foam::pow3(Cu[Kj+1])
              + alpha[2]*Foam::sqr(Cu[Kj+1])
              + alpha[1]*Cu[Kj+1]
              + alpha[0]
            )/6.0;

        const scalar Wmax =
            (
                alpha[3]*Foam::pow3(Cu[Kj])
              + alpha[2]*Foam::sqr(Cu[Kj])
              + alpha[1]*Cu[Kj]
              + alpha[0]
            )/6.0;

        if (Foam::mag(Wmin - Vf)/V < tol_)
        {
            // We are exactly on the upper bound

            return -Cu[Kj+1];
        }
        else if (Foam::mag(Wmax - Vf)/V < tol_)
        {
            // We are exactly on the lower bound

            return -Cu[Kj];
        }
        else if (Wmin <= Vf && Vf <= Wmax)
        {
            // Target falls within bracket

            const scalar z = LVE::exactCubicSolver
            (
                Vf,
                alpha,
                Cu[Kj],
                Cu[Kj+1],
                Wmin,
                Wmax,
                V
            );

            if (z != z)
                FatalErrorInFunction
                    << "LVE solver failed. C: " << z
                    << endl << abort(FatalError);

            return z;
        }
        else if (Vf > Wmax)
        {
            // Target volume is above bracket

            if (Kmin == Kj)
            {
                return -Cu[Kj];
            }
            else
            {
                // Adjust bracket
                Vmin = Wmax;
                Kmax = Kj;
            }
        }
        else
        {
            // Target volume is below bracket

            if (Kmax == Kj+1)
            {
                return -Cu[Kj+1];
            }
            else
            {
                // Adjust bracket
                Vmax = Wmin;
                Kmin = Kj+1;
            }
        }
    }

    FatalErrorInFunction
        <<  "LVE solver failed." << endl << abort(FatalError);

    return Cu[0];
}

inline scalar LVE::cLVE
(
    const vertexVector& v,
    const vector& n,
    const scalar f
) const
{
    // Consistently compute the volume of a general hexahedron and then call the
    // CIBRAVE algorithm

    vectorList triAreas(12);
    vectorList triNormals(12);
    scalarList C(12);
    labelList m(12);

    for (int j = 0; j < 12; j++)
    {
        triAreas[j] =
            (
                v[faceTriangleDecompCC[j][1]]
              - v[faceTriangleDecompCC[j][0]]
            )
          ^ (
                v[faceTriangleDecompCC[j][2]]
              - v[faceTriangleDecompCC[j][0]]
            );

        triNormals[j] = triAreas[j]/Foam::mag(triAreas[j]);

        C[j] = -(triNormals[j] & v[faceTriangleDecompCC[j][0]]);

        m[j] =
            Foam::mag(triNormals[j].y()) > Foam::mag(triNormals[j].x())
          ? Foam::mag(triNormals[j].z()) > Foam::mag(triNormals[j].y()) ? 2 : 1
          : Foam::mag(triNormals[j].z()) > Foam::mag(triNormals[j].x()) ? 2 : 0;
    }

    // Surface integration to calculate the volume

    scalar V = 0.0;
    for (int j = 0; j < 12; j++)
        V += -C[j]*triAreas[j][m[j]]/triNormals[j][m[j]]/6.0;

    return cLVE(v,n,f,V,triAreas,triNormals,C,m);
}

inline scalar LVE::pLVE
(
    const vertexVector& vertices,
    const vector& no,
    const scalar fo
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
        C = CMin + alpha*Foam::mag(CMax - CMin);
    }
    else
    {
        C = CMax - alpha*Foam::mag(CMax - CMin);
    }

    return C*S - (no & vertices.lba());
}

}

}

}
