#include "LVE.H"
#include "vof.H"
#include "truncatedHex.H"
#include "truncatedPiped.H"
#include "SortableList.H"
#include "constants.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

const scalar LVE::tol_ = 1e-10;
const label LVE::maxIter_ = 100;

scalar LVE::solveBracketedLVE
(
    const vertexVector& v,
    const scalar& V,
    const vector& n,
    const scalar& f,
    scalar CMin,
    scalar CMax,
    scalar VMin,
    scalar VMax
) const
{
    // Here we should exploit the analytical solution method of Lopez &
    // Hernandez (JCP, 2008). However, for now we've implemented a simple
    // bisection method assuming linear behavior between CMin and CMax.

    label iter = 0;

    const scalar Vf = f*V;

    scalar Vi = -V;
    scalar C = (CMin+CMax)/2.0;
    scalar pos = 0.5;

    while (Foam::mag((Vf-Vi)/V) > tol_ && iter < maxIter_)
    {
        // Use bisection of we are close to the edges of the [VMin,VMax]
        // interval. Use a linear interpolation otherwise.

        if (pos < 1e-3 || pos > 0.999)
        {
            C = (CMin+CMax)/2.0;
        }
        else
        {
            C = CMin + (Vf-VMin)/(VMax-VMin) * (CMax-CMin);
        }

        Vi = truncatedHex(v,n,C).volume();
        pos = (Vi - VMin)/(VMax-VMin+1e-50);

        #ifdef FULLDEBUG

        if (Vi - VMin < -1e-12 || VMax - Vi < -1e-12)
        {
            FatalErrorInFunction
                << "Bracketed LVE problem solver failed because the "
                << "truncated hex volume is non-monotonic." << endl
                << abort(FatalError);
        }

        #endif

        // Adjust search boundaries

        if (Vi < Vf)
        {
            CMin = C;
            VMin = Vi;
        }
        else
        {
            CMax = C;
            VMax = Vi;
        }

        iter++;
    }

    #ifdef FULLDEBUG

    if (iter == maxIter_)
    {
        WarningInFunction
            << "Bracketed LVE problem solver did not converge" << endl;
    }

    #endif

    return C;
}

LVE::LVE(const vof& vf)
:
    vf_(vf),
    fvMsh_(vf.fvMsh()),
    rectilinear_(fvMsh_.msh()[0].rectilinear() == unitXYZ)
{}

LVE::~LVE()
{}

scalar LVE::operator()(const labelVector& ijk, const vector& n) const
{
    const scalar f(vf_.alpha()[0][0](ijk));

    const vertexVector vertices
    (
        fvMsh_.template metrics<colocated>().vertexCenters()[0][0](ijk)
    );

    const vector m(n/Foam::mag(n));

    return rectilinear_ ? pLVE(vertices,m,f) : cLVE(vertices,m,f);
}

scalar LVE::aLVE
(
    const vertexVector& v,
    const vector& n,
    const scalar& f
) const
{
    scalarList C(8);
    scalarList V(8);

    for (int i = 0; i < 8; i++)
    {
        C[i] = -(n & v[i]);
        V[i] = truncatedHex(v,n,C[i]).volume();
    }

    labelList order;
    uniqueOrder(V, order);

    const scalar Vc = V[order[order.size()-1]];

    // Early detection at the vertices

    forAll(order, i)
    {
        if (Foam::mag(f-V[order[i]]/Vc) < tol_)
        {
            return C[order[i]];
        }
    }

    const scalar Vf = Vc*f;

    // Determine bracket

    scalar CMin = C[order[0]];
    scalar CMax = C[order[order.size()-1]];
    scalar VMin = V[order[0]];
    scalar VMax = V[order[order.size()-1]];

    for (int i = 1; i < order.size(); i++)
    {
        if (V[order[i-1]] < Vf && Vf < V[order[i]])
        {
            CMin = C[order[i-1]];
            CMax = C[order[i]];

            VMin = V[order[i-1]];
            VMax = V[order[i]];

            break;
        }
    }

    // Solve

    return solveBracketedLVE(v,Vc,n,f,CMin,CMax,VMin,VMax);
}

scalar LVE::cLVE
(
    const vertexVector& v,
    const vector& n,
    const scalar& f
) const
{

    /*

    +++++++++++++++++++++++++++++++++++++++++++++++

    CIBRAVE algorithm of J. López, J. Hernández,
    P. Gómez and F. Faura (2016) adapted for
    non-convex polyhedrons following . López,
    J. Hernández, P. Gómez and F. Faura (2019)

    +++++++++++++++++++++++++++++++++++++++++++++++

    */

    scalar Vmin = 0, Vmax;
    scalar Vkmin, Vkmax;
    int Kmin = 0, Kmax = 7, Kf = 0;
    int maxNumberOfFaces = 16;

    scalarList Ajs(12);
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
    int NIPV0[12] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    int IPV0[12][5] = {
        {2, 0, 6, 2, 0},
        {6, 0, 4, 6, 0},
        {5, 1, 3, 5, 1},
        {5, 3, 7, 5, 3},
        {5, 0, 1, 5, 0},
        {4, 0, 5, 4, 0},
        {6, 3, 2, 6, 3},
        {7, 3, 6, 7, 3},
        {3, 0, 2, 3, 0},
        {1, 0, 3, 1, 0},
        {6, 4, 5, 6, 4},
        {6, 5, 7, 6, 5}
    };

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

    vector aux;

    for (int j = 0; j < 12; j++)
    {
        aux = (v[IPV0[j][2]] - v[IPV0[j][1]])
            ^ (v[IPV0[j][3]] - v[IPV0[j][1]]);
        ns[j] = aux / Foam::mag(aux);
        Cjs[j] =  - ns[j] & v[IPV0[j][1]];
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

    for (int j = 0; j < 12; j++)
    {
        hjs[j] = int((NIPV0[j]-2)/2);
        epsilonjs[j] = NIPV0[j] % 2;

        Ajs[j] = Foam::zero();

        for (int i = 2; i <= (hjs[j] + 1); i++)
        {
            Ajs[j] += vectorialProductProjection(
                (v[IPV0[j][2*i - 1]] - v[IPV0[j][1]]),
                (v[IPV0[j][2*i]] - v[IPV0[j][2*i-2]]),
                maxNormalIndex[j]
                );
        }

        if (epsilonjs[j] == 1)
        {
            Ajs[j] += vectorialProductProjection(
                (v[IPV0[j][NIPV0[j]]] - v[IPV0[j][1]]),
                (v[IPV0[j][1]] - v[IPV0[j][NIPV0[j]-1]]),
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

    Cs[0] = SortedCs[0];
    Kmap[sortedCsIndex[0]] = 0;
    int Knum = 0;

    for (int i = 1; i < 8; i++)
    {
        if (Foam::mag(SortedCs[i-1] - SortedCs[i]) < 1e-12)
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

        ////////////////////////////////////

        // Find Bracket by interpolation

        ////////////////////////////////////

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

        ////////////////////////////////////

        // Arrange New Face

        ////////////////////////////////////

        for (int j = 0; j < 12; j++)
        {
            for (int i = 1; i <= 3; i++)
            {
                if (IA[IPV0[j][i]] == 1)
                    NIPV1[j]++;
            }

            if (NIPV1[j] > 0)
            {
                int index = 0;
                for (int i = 1; i <= NIPV0[j]; i++)
                {
                    int ip1 = IPV0[j][i];
                    int ip2 = IPV0[j][i+1];

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

        ////////////////////////////////////

        // Compute the picewise analytical function

        ////////////////////////////////////

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
                    Kjs[j] += vectorialProductProjection(
                        (x0[IPV1[j][2*i - 1]] - x0[IPV1[j][1]]),
                        (x0[IPV1[j][2*i]] - x0[IPV1[j][2*i-2]]),
                        maxNormalIndex[j]
                        );
                    Ljs[j] += vectorialProductProjection(
                        (x0[IPV1[j][2*i - 1]] - x0[IPV1[j][1]]),
                        (E[IPV1[j][2*i]] - E[IPV1[j][2*i-2]]),
                        maxNormalIndex[j]
                        )
                        + vectorialProductProjection(
                        (E[IPV1[j][2*i - 1]] - E[IPV1[j][1]]),
                        (x0[IPV1[j][2*i]] - x0[IPV1[j][2*i-2]]),
                        maxNormalIndex[j]
                        );
                    Mjs[j] += vectorialProductProjection(
                        (E[IPV1[j][2*i - 1]] - E[IPV1[j][1]]),
                        (E[IPV1[j][2*i]] - E[IPV1[j][2*i-2]]),
                        maxNormalIndex[j]
                        );
                }

                if (epsilonjs[j] == 1)
                {
                    Kjs[j] += vectorialProductProjection(
                        (x0[IPV1[j][NIPV1[j]]] - x0[IPV1[j][1]]),
                        (x0[IPV1[j][1]] - x0[IPV1[j][NIPV1[j]-1]]),
                        maxNormalIndex[j]
                        );
                    Ljs[j] += vectorialProductProjection(
                        (x0[IPV1[j][NIPV1[j]]] - x0[IPV1[j][1]]),
                        (E[IPV1[j][1]] - E[IPV1[j][NIPV1[j]-1]]),
                        maxNormalIndex[j]
                        )
                        + vectorialProductProjection(
                        (E[IPV1[j][NIPV1[j]]] - E[IPV1[j][1]]),
                        (x0[IPV1[j][1]] - x0[IPV1[j][NIPV1[j]-1]]),
                        maxNormalIndex[j]
                        );
                    Mjs[j] += vectorialProductProjection(
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

        ////////////////////////////////////

        // Check bracket and solve inverse problem

        ///////////////////////////////////

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

scalar LVE::pLVE
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

scalar LVE::exactCubicSolver
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
    const scalar Cr = Foam::mag(Cmax - Cmin);
    if (Foam::mag(alpha3) > 0)
    {
        if (Foam::min(
                Foam::min(
                    Foam::mag(alpha3/alpha0),
                    Foam::mag(alpha3/alpha1)),
                Foam::mag(alpha3/alpha2)
            ) < tol_)
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
                const scalar phi2 = phi1 - ((2*Foam::constant::mathematical::pi)/3);
                const scalar phi3 = phi1 + ((2*Foam::constant::mathematical::pi)/3);

                const scalar z1 = 2 * Foam::sqrt(-q) * Foam::cos(phi1) - a2 /3;
                const scalar z2 = 2 * Foam::sqrt(-q) * Foam::cos(phi2) - a2 /3;
                const scalar z3 = 2 * Foam::sqrt(-q) * Foam::cos(phi3) - a2 /3;

                if ((Cmin <= z1 + (tol_ * Cr)) && (z1 <= Cmax + (tol_ * Cr)))
                {
                    return -z1;
                }
                else if ((Cmin <= z2 + (tol_ * Cr)) && (z2 <= Cmax + (tol_ * Cr)))
                {
                    return -z2;
                }
                else if ((Cmin <= z3 + (tol_ * Cr)) && (z3 <= Cmax + (tol_ * Cr)))
                {
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
                    return -z1;
                }
                else if ((Cmin <= z2 + (tol_ * Cr)) && (z2 <= Cmax + (tol_ * Cr)))
                {
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


scalar LVE::newtonCubicSolver
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

    scalar Cn = Cmin + (Cmax - Cmin) * (Vf - Vkmax)/(Vkmin - Vkmax);
    scalar fn, dfn;

    for (int iter = 0; iter < maxIter_; iter++)
    {
        fn = (1/6.0) * ((alpha3 * Foam::pow3(Cn))
            + (alpha2 * Foam::sqr(Cn)) + (alpha1 * Cn) + alpha0) - Vf;

        if (Foam::mag(fn/Vtotal) < tol_)
        {
            return - Cn;
        }

        dfn = (1/6.0) * ((3.0 * alpha3 * Foam::sqr(Cn))
            + (2.0 * alpha2 * Cn) + (alpha1));

        Cn = Cn - fn / dfn;
    }

    Info << "Warning!! LVE solver Newton method didn't converge! Last error:  " << fn/Vtotal << endl;
    return -Cn;

}

}

}

}


