#include "LSFIR.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"
#include "LVE.H"
#include "SortableList.H"
#include "constants.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(LSFIR, 0);
addToRunTimeSelectionTable(normalScheme, LSFIR, dictionary);

LSFIR::LSFIR(const vof& vf, const dictionary& dict)
:
    normalScheme(vf, dict),
    threshold_(dict.lookupOrDefault<scalar>("threshold", 1e-3))
{}

LSFIR::LSFIR(const LSFIR& s)
:
    normalScheme(s),
    threshold_(s.threshold_)
{}

LSFIR::~LSFIR()
{}

tmp<colocatedVectorField> LSFIR::operator()()
{


    tmp<colocatedVectorField> tn
    (
        new colocatedVectorField(
            "normal",
            vf_.fvMsh()
        )
    );

    colocatedVectorDirection& n = tn.ref()[0][0];

    colocatedScalarDirection alpha = vf_.alpha()[0][0];
    const colocatedVertexVectorDirection& v =
        vf_.fvMsh().template
            metrics<colocated>().vertexCenters()[0][0];
    const meshDirection<vector,colocated>& centers =
        vf_.fvMsh().template
            metrics<colocated>().cellCenters()[0][0];

    tmp<colocatedVectorField> xn
    (
        new colocatedVectorField(
            "surface_centers",
            vf_.fvMsh()
        )
    );

    colocatedVectorDirection& xgi = xn.ref()[0][0];

    /*

    Initialise the field using LSGIR

    */

    forAllCells(n, i, j, k)
    {
        if ((alpha(i,j,k) > 1e-12) && (alpha(i,j,k) < 1 - 1e-12))
        {

            double Aaux[26][3];
            double baux[26];
            int index;
            scalarList weights(26);
            tensor A = Zero;
            vector b = Zero;

            for (int aux1 = 0; aux1 < 3; aux1++)
            {
                for (int aux2 = 0; aux2 < 3; aux2++)
                {
                    for (int aux3 = 0; aux3 < 3; aux3++)
                    {
                        if ((aux1 != 1) || (aux2 != 1) || (aux3 != 1))
                        {
                            index = aux1 + 3 * aux2 + 9 * aux3;
                            if (index > 13)
                                index--;
                            weights[index] = (1.0) / Foam::sqr(Foam::mag(centers(i,j,k) - centers(i+aux1-1,j+aux2-1,k+aux3-1)));
                            baux[index] = weights[index] * (alpha(i+aux1-1,j+aux2-1,k+aux3-1) - alpha(i,j,k));
                            Aaux[index][0] = weights[index] * (centers(i+aux1-1,j+aux2-1,k+aux3-1)[0] - centers(i,j,k)[0]);
                            Aaux[index][1] = weights[index] * (centers(i+aux1-1,j+aux2-1,k+aux3-1)[1] - centers(i,j,k)[1]);
                            Aaux[index][2] = weights[index] * (centers(i+aux1-1,j+aux2-1,k+aux3-1)[2] - centers(i,j,k)[2]);
                        }
                    }
                }
            }

            for (int aux3 = 0; aux3 < 26; aux3++)
            {
                for (int aux1 = 0; aux1 < 3; aux1++)
                {
                    for (int aux2 = 0; aux2 < 3; aux2++)
                    {
                        A[3 * aux1 + aux2] += Aaux[aux3][aux1] * Aaux[aux3][aux2];
                    }
                    b[aux1] += Aaux[aux3][aux1] * baux[aux3];
                }
            }

            n(i,j,k) = A.inv() & b;
            const scalar S = Foam::mag(n(i,j,k));
            n(i,j,k) /= S;

        }
        else
        {
            n(i,j,k) = Zero;
        }
    }

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
    int originalNf = 12;
    int maxIterNumber = 4;
    bool updatedNormals = false;

    for (int iter = 0; iter < maxIterNumber; iter ++)
    {
        if (updatedNormals)
        {
            break;
        }

        updatedNormals = true;

        forAllCells(n, i, j, k)
        {
            xgi(i,j,k) = Foam::zero();

            if ((alpha(i,j,k) > 1e-12) && (alpha(i,j,k) < (1 - 1e-12)))
            {

                const scalar S = Foam::mag(n(i,j,k));

                if (S < threshold_)
                {
                    n(i,j,k)[0] = 1;
                }

                const scalar C = vf_.lve()(i,j,k,n(i,j,k));
                double originalCs[8];

                for (int il = 0; il < 8; il++)
                {
                    originalCs[il] = (n(i,j,k) & v(i,j,k)[il]);
                }

                double SortedCs[8] = {0};
                int sortedCsIndex[8] = {0};

                int len = 1;
                int place = 0;

                SortedCs[0] = originalCs[0];
                sortedCsIndex[0] = 0;

                for (int il = 1; il < 8; il++)
                {
                    place = len;
                    for (int jl = 0; jl < len; jl++)
                    {
                        if (originalCs[il] < SortedCs[jl])
                        {
                            place = jl;
                            break;
                        }
                    }

                    for (int jl = len; jl > place; jl--)
                    {
                    SortedCs[jl] = SortedCs[jl-1];
                    sortedCsIndex[jl] = sortedCsIndex[jl-1];
                    }

                    SortedCs[place] = originalCs[il];
                    sortedCsIndex[place] = il;
                    len++;

                }

                int Kmap[8];
                scalarList Cs(8);
                int IA[8];
                int Kf = 0;

                Cs[0] = SortedCs[0];
                Kmap[sortedCsIndex[0]] = 0;
                int Knum = 0;

                for (int il = 1; il < 8; il++)
                {
                    if (Foam::mag(SortedCs[il-1] - SortedCs[il]) < 1e-12)
                    {
                        Kmap[sortedCsIndex[il]] = Knum;
                    }
                    else
                    {
                        Knum++;
                        Cs[Knum] = SortedCs[il];
                        Kmap[sortedCsIndex[il]] = Knum;
                    }
                }

                for (int il = 1; il <= Knum; il++)
                {
                    if (Cs[il] > -C)
                    {
                        Kf = il - 1;
                        break;
                    }

                    if (il == Knum)
                    {
                        Kf = il - 1;
                    }
                }


                for (int il = 0; il < 8; il++)
                {
                    if (Kmap[il] <= Kf)
                        IA[il] = 0;
                    else
                        IA[il] = 1;
                }

                int NIPV1[16] = {0};
                int IPV1[16][14];

                int insertedVertex[12][6];
                int totalInsertedVertex = 0;
                int Nf = originalNf;
                int Ip = 8;

                for (int jt = 0; jt < originalNf; jt++)
                {
                    for (int it = 1; it <= NIPV0[jt]; it++)
                    {
                        if (IA[IPV0[jt][it]] == 1)
                            NIPV1[jt]++;
                    }

                    if (NIPV1[jt] > 0)
                    {
                        int index = 0;
                        for (int it = 1; it <= NIPV0[jt]; it++)
                        {
                            int ip1 = IPV0[jt][it];
                            int ip2 = IPV0[jt][it+1];

                            if (IA[ip1] == 1)
                            {
                                index++;
                                IPV1[jt][index] = ip1;
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
                                        IPV1[jt][index] = insertedVertex[aux][0];
                                        if (IA[ip2] == 1)
                                        {
                                            insertedVertex[aux][3] = jt;
                                            insertedVertex[aux][4] = index;
                                        }
                                        break;
                                    }
                                }

                                if (flag)
                                {
                                    insertedVertex[totalInsertedVertex][0] = Ip;
                                    IPV1[jt][index] = Ip;

                                    if (IA[ip2] == 1)
                                    {
                                        insertedVertex[totalInsertedVertex][3] = jt;
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

                        NIPV1[jt] = index;
                        IPV1[jt][0] = IPV1[jt][NIPV1[jt]];
                        IPV1[jt][NIPV1[jt]+1] = IPV1[jt][1];
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

                if (Nf - originalNf == 1)
                {
                    vector Xin;
                    vector e;
                    vectorList x0(totalInsertedVertex);

                    for (int it = 0; it < totalInsertedVertex; it++)
                    {
                        e = v(i,j,k)[insertedVertex[it][2]] - v(i,j,k)[insertedVertex[it][1]];
                        e /= Foam::mag(e);
                        Xin = v(i,j,k)[insertedVertex[it][1]];
                        x0[insertedVertex[it][0] - 8] = Xin + (-C - (n(i,j,k) & Xin)) * e * (1 / (n(i,j,k) & e));
                    }

                    scalar auxMax;
                    int maxNormalIndex;

                    auxMax = Foam::mag(n(i,j,k)[0]);
                    maxNormalIndex = 0;

                    for (int it = 1; it < 3; it++)
                    {
                        if (Foam::mag(n(i,j,k)[it]) > auxMax)
                        {
                            auxMax = Foam::mag(n(i,j,k)[it]);
                            maxNormalIndex = it;
                        }
                    }

                    scalar TotalArea = 0;
                    vector Centroid;
                    scalar Area;

                    for (int it = 2; it < NIPV1[originalNf]; it++)
                    {
                        Centroid = x0[IPV1[originalNf][1]-8];
                        Centroid += x0[IPV1[originalNf][it]-8];
                        Centroid += x0[IPV1[originalNf][it+1]-8];
                        Centroid /= 3;

                        Area = vectorialProductProjection(
                            x0[IPV1[originalNf][it]-8] - x0[IPV1[originalNf][1]-8],
                            x0[IPV1[originalNf][it+1]-8] - x0[IPV1[originalNf][1]-8],
                            maxNormalIndex);
                        Area /= n(i,j,k)[maxNormalIndex];

                        TotalArea += Area;
                        xgi(i,j,k) += Area * Centroid;
                    }

                    xgi(i,j,k) /= TotalArea;

                    if (xgi(i,j,k)[1] != xgi(i,j,k)[1])
                    {
                    FatalErrorInFunction
                        <<  "LSFIR: geometric center computation failed." << endl << abort(FatalError);
                    }
                }
                else
                {
                    int numberOfPoints = 0;
                    vector e;
                    vector Xin;

                    for (int aux1 = 0; aux1 < originalNf; aux1++)
                    {
                        for (int aux2 = 1; aux2 <= NIPV0[aux1]; aux2++)
                        {
                            if (((IA[IPV0[aux1][aux2]] == 0) && (IA[IPV0[aux1][aux2 + 1]] == 1)) ||
                                ((IA[IPV0[aux1][aux2]] == 1) && (IA[IPV0[aux1][aux2 + 1]] == 0)))
                            {
                                if (IA[IPV0[aux1][aux2]] == 1)
                                {
                                    e = v(i,j,k)[IPV0[aux1][aux2 + 1]] - v(i,j,k)[IPV0[aux1][aux2]];
                                    Xin = v(i,j,k)[IPV0[aux1][aux2]];
                                }
                                else
                                {
                                    e = v(i,j,k)[IPV0[aux1][aux2]] - v(i,j,k)[IPV0[aux1][aux2 + 1]];
                                    Xin = v(i,j,k)[IPV0[aux1][aux2 + 1]];
                                }

                                e /= Foam::mag(e);
                                xgi(i,j,k) += Xin + (-C - (n(i,j,k) & Xin)) * e * (1 / (n(i,j,k) & e));
                                numberOfPoints++;
                            }
                        }
                    }
                    xgi(i,j,k) /= numberOfPoints;
                }

                if (xgi(i,j,k)[1] != xgi(i,j,k)[1])
                {
                    FatalErrorInFunction
                        <<  "LSFIR: center of vertex computation failed." << endl << abort(FatalError);
                }
            }
        }

        xn.ref()[0].correctBoundaryConditions(true);

        forAllCells(n, i, j, k)
        {

            if ((alpha(i,j,k) > 1e-12) && (alpha(i,j,k) < 1 - 1e-12))
            {

                int d1 = 0;
                int d2 = 1;
                int d3 = 2;
                scalar value;
                double A[2][2] = {0};
                double b[2] = {0};
                scalar wi = 0;

                value = Foam::mag(n(i,j,k)[0]);

                if (Foam::mag(n(i,j,k)[1]) > value)
                {
                    d1 = 1;
                    d2 = 0;
                    value = Foam::mag(n(i,j,k)[1]);
                }

                if (Foam::mag(n(i,j,k)[2]) > value)
                {
                    d3 = d1;
                    d1 = 2;
                    value = Foam::mag(n(i,j,k)[2]);
                }

                for (int aux1 = 0; aux1 < 3; aux1++)
                {
                    for (int aux2 = 0; aux2 < 3; aux2++)
                    {
                        for (int aux3 = 0; aux3 < 3; aux3++)
                        {
                            if (
                                (alpha(i+aux1-1,j+aux2-1,k+aux3-1) > 1e-12) &&
                                (alpha(i+aux1-1,j+aux2-1,k+aux3-1) < 1 - 1e-12) &&
                                ((aux1 != 1) || (aux2 != 1) || (aux3 != 1))
                                )
                            {
                                scalar angle = (n(i+aux1-1,j+aux2-1,k+aux3-1) & n(i,j,k))
                                                / (Foam::mag(n(i+aux1-1,j+aux2-1,k+aux3-1)) * Foam::mag(n(i,j,k)));

                                if (angle > Foam::cos(0.25 * Foam::constant::mathematical::pi))
                                {
                                    wi = (1.0)/Foam::pow(Foam::mag(xgi(i+aux1-1,j+aux2-1,k+aux3-1)-xgi(i,j,k)),2.5);
                                    A[0][0] += wi*Foam::sqr(xgi(i+aux1-1,j+aux2-1,k+aux3-1)[d2] - xgi(i,j,k)[d2]);
                                    A[0][1] += wi*(xgi(i+aux1-1,j+aux2-1,k+aux3-1)[d2] - xgi(i,j,k)[d2])
                                        *(xgi(i+aux1-1,j+aux2-1,k+aux3-1)[d3] - xgi(i,j,k)[d3]);
                                    A[1][1] += wi*Foam::sqr(xgi(i+aux1-1,j+aux2-1,k+aux3-1)[d3] - xgi(i,j,k)[d3]);
                                    b[0] += wi*(xgi(i+aux1-1,j+aux2-1,k+aux3-1)[d2] - xgi(i,j,k)[d2])
                                        *(xgi(i+aux1-1,j+aux2-1,k+aux3-1)[d1] - xgi(i,j,k)[d1]);
                                    b[1] += wi*(xgi(i+aux1-1,j+aux2-1,k+aux3-1)[d3] - xgi(i,j,k)[d3])
                                        *(xgi(i+aux1-1,j+aux2-1,k+aux3-1)[d1] - xgi(i,j,k)[d1]);
                                }
                            }
                        }
                    }
                }

                A[1][0] = A[0][1];

                vector nc;

                if ((A[0][0] * A[1][1] - A[0][1] * A[1][0]) != 0)
                {
                    nc[d1] = n(i,j,k)[d1] > 0 ? 1 : -1;
                    nc[d2] = (A[1][1] * b[0] - A[0][1] * b[1]) / (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
                    nc[d3] = (- A[1][0] * b[0] + A[0][0] * b[1]) / (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
                    nc /= Foam::mag(nc);
                    scalar angle = (nc & n(i,j,k))
                                    / (Foam::mag(nc) * Foam::mag(n(i,j,k)));
                    if (angle > Foam::cos((1.0 / 6.0) * Foam::constant::mathematical::pi / double(iter)))
                    {
                        n(i,j,k) = nc;
                        updatedNormals = false;
                    }
                }

                if (n(i,j,k)[1] != n(i,j,k)[1])
                {
                    FatalErrorInFunction
                        <<  "LSFIR: normal computation failed." << endl << abort(FatalError);
                }
            }
        }
    }

    return tn;
}

}

}

}
