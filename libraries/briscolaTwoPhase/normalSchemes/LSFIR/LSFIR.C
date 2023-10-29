#include "LSFIR.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"
#include "LVE.H"
#include "SortableList.H"
#include "constants.H"
#include "Time.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(LSFIR, 0);
addToRunTimeSelectionTable(normalScheme, LSFIR, dictionary);

void LSFIR::createBoundaryTypes()
{
    const faceLabel& faceType = fvMsh_.msh().facePatchType();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                bool aux2 = true;

                if (i != 1)
                {
                    aux2 =
                        (i == 0)
                      ? (aux2 && faceType.left() > 0)
                      : (aux2 && faceType.right() > 0);
                }

                if (j != 1)
                {
                    aux2 =
                        (j == 0)
                      ? (aux2 && faceType.bottom() > 0)
                      : (aux2 && faceType.top() > 0);
                }

                if (k != 1)
                {
                    aux2 =
                        (k == 0)
                      ? (aux2 && faceType.aft() > 0)
                      : (aux2 && faceType.fore() > 0);
                }

                boundaryTypeLSFIR_[i][j][k] = aux2;
            }
        }
    }
}

void LSFIR::createHexagonDescription()
{
    // Create the table to describe the hexahedron for planar faces or
    // non-planar faces hexahedrons

    if (planarFaces_)
    {
        int NIPV0[12] = {4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0};

        for (int i = 0; i < 6; i++)
        {
            NIPV0_[i] = NIPV0[i];
            for (int j = 0; j < 4; j++)
            {
                IPV0_[i][j + 1] = vertexNumsInFaceCC[i][j];
            }

            IPV0_[i][0] = IPV0_[i][NIPV0[i]];
            IPV0_[i][NIPV0[i] + 1] = IPV0_[i][1];
        }

        originalNf_ = 6;
    }
    else
    {
        int NIPV0[12] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

        for (int i = 0; i < 12; i++)
        {
            NIPV0_[i] = NIPV0[i];
            for (int j = 0; j < 3; j++)
            {
                IPV0_[i][j + 1] = faceTriangleDecompCC[i][j];
            }

            IPV0_[i][0] = IPV0_[i][NIPV0[i]];
            IPV0_[i][NIPV0[i] + 1] = IPV0_[i][1];
        }

        originalNf_ = 12;
    }
}

LSFIR::LSFIR
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
:
    LSGIR(fvMsh, dict, alpha),
    planarFaces_(dict.lookupOrDefault<bool>("planarFaces", true)),
    centerAveraging_(dict.lookupOrDefault<bool>("centerAveraging", false)),
    lve_(fvMsh.msh()[0].rectilinear() == unitXYZ)
{
    createBoundaryTypes();
    createHexagonDescription();
}

LSFIR::LSFIR(const LSFIR& s)
:
    LSGIR(s),
    planarFaces_(s.planarFaces_),
    centerAveraging_(s.centerAveraging_),
    lve_(s.fvMsh().msh()[0].rectilinear() == unitXYZ)
{
    createBoundaryTypes();
    createHexagonDescription();
}

LSFIR::~LSFIR()
{}

void LSFIR::correct()
{
    colocatedVectorField& n = *this;

    // Reconstruct the interface normal using the version of LSFIR presented in
    // Lopez (2022).

    const colocatedVertexVectorField& v =
        fvMsh_.template metrics<colocated>().vertexCenters();
    const colocatedScalarField& cv =
        fvMsh_.template metrics<colocated>().cellVolumes();

    colocatedVectorField xgi("surface_centers", fvMsh_);

    const scalar angleTol = Foam::cos(1e-3);

    const scalar validAngleTol =
        Foam::cos(0.25 * Foam::constant::mathematical::pi);

    // Initialise the normal using LSGIR and store that one

    LSGIR::correct();
    colocatedVectorField nNew(n);

    int maxIterNumber = 4;
    int updatedNormals = 1;

    for (int iter = 0; iter < maxIterNumber; iter ++)
    {
        // Check the stopping criterion

        reduce(updatedNormals, maxOp<int>());

        if (updatedNormals == 0)
        {
            break;
        }

        scalar maxAngleTol = Foam::cos
        (
            0.166 * Foam::constant::mathematical::pi
          / scalar(iter + 1)
        );

        updatedNormals = 0;

        // Compute the centers of the interface for all interfacial cells

        forAllCells(n, i, j, k)
        {
            xgi(i,j,k) = Zero;

            if
            (
                alpha_(i,j,k) > vof::threshold
             && alpha_(i,j,k) < (1.0 - vof::threshold)
            )
            {

                // Prevent computing interface for zero normal

                const scalar S = Foam::mag(n(i,j,k));

                if (S < 1e-12)
                {
                    n(i,j,k)[0] = 1;
                }

                // Compute which vertex are inside the liquid for the truncated
                // cell

                scalar C =
                    lve_(alpha_(i,j,k),v(i,j,k),cv(i,j,k),n(i,j,k));

                double originalCs[8];

                for (int aux1 = 0; aux1 < 8; aux1++)
                {
                    originalCs[aux1] = (n(i,j,k) & v(i,j,k)[aux1]);
                }

                double SortedCs[8] = {0};
                int sortedCsIndex[8] = {0};

                int len = 1;
                int place = 0;

                SortedCs[0] = originalCs[0];
                sortedCsIndex[0] = 0;

                for (int aux1 = 1; aux1 < 8; aux1++)
                {
                    place = len;
                    for (int aux2 = 0; aux2 < len; aux2++)
                    {
                        if (originalCs[aux1] < SortedCs[aux2])
                        {
                            place = aux2;
                            break;
                        }
                    }

                    for (int aux2 = len; aux2 > place; aux2--)
                    {
                    SortedCs[aux2] = SortedCs[aux2-1];
                    sortedCsIndex[aux2] = sortedCsIndex[aux2-1];
                    }

                    SortedCs[place] = originalCs[aux1];
                    sortedCsIndex[place] = aux1;
                    len++;

                }

                int Kmap[8];
                scalarList Cs(8);
                int IA[8];
                int Kf = 0;

                Cs[0] = SortedCs[0];
                Kmap[sortedCsIndex[0]] = 0;
                int Knum = 0;

                for (int aux1 = 1; aux1 < 8; aux1++)
                {
                    if
                    (
                        Foam::mag(SortedCs[aux1-1] - SortedCs[aux1])
                      / Foam::mag(SortedCs[0] - SortedCs[7])
                      < 1e-8
                    )
                    {
                        Kmap[sortedCsIndex[aux1]] = Knum;
                    }
                    else
                    {
                        Knum++;
                        Cs[Knum] = SortedCs[aux1];
                        Kmap[sortedCsIndex[aux1]] = Knum;
                    }
                }

                if (Foam::mag(Cs[0] + C)/(Cs[Knum] - Cs[0]) < 1e-8)
                {
                    C += - (Cs[Knum] - Cs[0]) * 1e-8;
                }
                else if (Foam::mag(Cs[Knum] + C)/(Cs[Knum] - Cs[0]) < 1e-8)
                {
                    C += (Cs[Knum] - Cs[0]) * 1e-8;
                }

                for (int aux1 = 1; aux1 <= Knum; aux1++)
                {
                    if (Cs[aux1] > -C)
                    {
                        Kf = aux1 - 1;
                        break;
                    }

                    if (aux1 == Knum)
                    {
                        Kf = aux1 - 1;
                    }
                }

                for (int aux1 = 0; aux1 < 8; aux1++)
                {
                    if (Kmap[aux1] <= Kf)
                        IA[aux1] = 0;
                    else
                        IA[aux1] = 1;
                }

                // If centerAveraging is true compute the center as the average
                // of the vertex, if not compute it as a center of mass

                if (!centerAveraging_)
                {
                    // Find the truncated interface using Lopez (2020) algorithm

                    int NIPV1[16] = {0};
                    int IPV1[16][14];

                    int insertedVertex[12][6];
                    int totalInsertedVertex = 0;
                    int Nf = originalNf_;
                    int Ip = 8;

                    for (int jt = 0; jt < originalNf_; jt++)
                    {
                        for (int it = 1; it <= NIPV0_[jt]; it++)
                        {
                            if (IA[IPV0_[jt][it]] == 1)
                                NIPV1[jt]++;
                        }

                        if (NIPV1[jt] > 0)
                        {
                            int index = 0;
                            for (int it = 1; it <= NIPV0_[jt]; it++)
                            {
                                int ip1 = IPV0_[jt][it];
                                int ip2 = IPV0_[jt][it+1];

                                if (IA[ip1] == 1)
                                {
                                    index++;
                                    IPV1[jt][index] = ip1;
                                }
                                if (IA[ip1] != IA[ip2])
                                {
                                    index++;
                                    bool flag = true;

                                    for
                                    (
                                        int aux = 0;
                                        aux < totalInsertedVertex;
                                        aux ++
                                    )
                                    {
                                        if
                                        (
                                            (
                                                (insertedVertex[aux][1] == ip1)
                                             || (insertedVertex[aux][2] == ip1)
                                            )
                                         && (
                                                (insertedVertex[aux][1] == ip2)
                                             || (insertedVertex[aux][2] == ip2)
                                            )
                                        )
                                        {
                                            flag = false;
                                            IPV1[jt][index] =
                                                insertedVertex[aux][0];

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
                                        insertedVertex
                                            [totalInsertedVertex][0] = Ip;

                                        IPV1[jt][index] = Ip;

                                        if (IA[ip2] == 1)
                                        {
                                            insertedVertex
                                                [totalInsertedVertex][3] = jt;
                                            insertedVertex
                                                [totalInsertedVertex][4] =
                                                    index;
                                            insertedVertex
                                                [totalInsertedVertex][2] = ip1;
                                            insertedVertex
                                                [totalInsertedVertex][1] = ip2;
                                        }
                                        else
                                        {
                                            insertedVertex
                                                [totalInsertedVertex][1] = ip1;
                                            insertedVertex
                                                [totalInsertedVertex][2] = ip2;
                                        }

                                        insertedVertex
                                            [totalInsertedVertex][5] = 0;

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
                                for
                                (
                                    int auxip = 0;
                                    auxip < totalInsertedVertex;
                                    auxip++
                                )
                                {
                                    if (insertedVertex[auxip][5] == 0)
                                    {
                                        IPV1[Nf][index] =
                                            insertedVertex[auxip][0];
                                        insertedVertex[auxip][5] = 1;
                                        taggedVertex++;
                                        break;
                                    }
                                }
                            }

                        }

                        counter++;
                    }

                    // If there are more than one interface polygons
                    // compute the center as an average of vertex

                    if (Nf - originalNf_ == 1)
                    {
                        vector Xin;
                        vector e;
                        vectorList x0(totalInsertedVertex);

                        for (int it = 0; it < totalInsertedVertex; it++)
                        {
                            e = v(i,j,k)[insertedVertex[it][2]]
                                - v(i,j,k)[insertedVertex[it][1]];
                            e /= Foam::mag(e);
                            Xin = v(i,j,k)[insertedVertex[it][1]];
                            x0[insertedVertex[it][0] - 8] =
                                Xin
                              + (-C - (n(i,j,k) & Xin))
                              * e * (1.0 / (n(i,j,k) & e));
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

                        for (int it = 2; it < NIPV1[originalNf_]; it++)
                        {
                            Centroid = x0[IPV1[originalNf_][1]-8];
                            Centroid += x0[IPV1[originalNf_][it]-8];
                            Centroid += x0[IPV1[originalNf_][it+1]-8];
                            Centroid /= 3;

                            Area =
                                crossProductComponent
                                (
                                    x0[IPV1[originalNf_][it]-8]
                                  - x0[IPV1[originalNf_][1]-8],
                                    x0[IPV1[originalNf_][it+1]-8]
                                  - x0[IPV1[originalNf_][1]-8],
                                    maxNormalIndex
                                );

                            Area /= n(i,j,k)[maxNormalIndex];

                            TotalArea += Area;
                            xgi(i,j,k) += Area * Centroid;
                        }

                        xgi(i,j,k) /= TotalArea;

                        if (xgi(i,j,k)[1] != xgi(i,j,k)[1])
                        {
                            FatalErrorInFunction
                                << "LSFIR: geometric center computation failed."
                                << endl << abort(FatalError);
                        }
                    }
                    else
                    {
                        int numberOfPoints = 0;
                        vector e;
                        vector Xin;

                        for (int aux1 = 0; aux1 < originalNf_; aux1++)
                        {
                            for (int aux2 = 1; aux2 <= NIPV0_[aux1]; aux2++)
                            {
                                if
                                (
                                    IA[IPV0_[aux1][aux2]] == 0
                                 && IA[IPV0_[aux1][aux2 + 1]] == 1
                                )
                                {
                                    if (IA[IPV0_[aux1][aux2]] == 1)
                                    {
                                        e =
                                            v(i,j,k)[IPV0_[aux1][aux2 + 1]]
                                          - v(i,j,k)[IPV0_[aux1][aux2]];

                                        Xin = v(i,j,k)[IPV0_[aux1][aux2]];
                                    }
                                    else
                                    {
                                        e =
                                            v(i,j,k)[IPV0_[aux1][aux2]]
                                          - v(i,j,k)[IPV0_[aux1][aux2 + 1]];

                                        Xin = v(i,j,k)[IPV0_[aux1][aux2 + 1]];
                                    }

                                    e /= Foam::mag(e);

                                    xgi(i,j,k) += Xin
                                        + (-C - (n(i,j,k) & Xin))
                                        * e * (1 / (n(i,j,k) & e));

                                    numberOfPoints++;
                                }
                            }
                        }

                        xgi(i,j,k) /= numberOfPoints;

                        if (xgi(i,j,k)[1] != xgi(i,j,k)[1])
                        {
                            FatalErrorInFunction
                                << "LSFIR: center of vertex computation failed."
                                << endl << abort(FatalError);
                        }
                    }
                }
                else
                {
                    int numberOfPoints = 0;
                    vector e;
                    vector Xin;

                    for (int aux1 = 0; aux1 < originalNf_; aux1++)
                    {
                        for (int aux2 = 1; aux2 <= NIPV0_[aux1]; aux2++)
                        {
                            if
                            (
                                IA[IPV0_[aux1][aux2]] == 0
                             && IA[IPV0_[aux1][aux2 + 1]] == 1
                            )
                            {
                                if (IA[IPV0_[aux1][aux2]] == 1)
                                {
                                    e =
                                        v(i,j,k)[IPV0_[aux1][aux2 + 1]]
                                      - v(i,j,k)[IPV0_[aux1][aux2]];

                                    Xin = v(i,j,k)[IPV0_[aux1][aux2]];
                                }
                                else
                                {
                                    e =
                                        v(i,j,k)[IPV0_[aux1][aux2]]
                                      - v(i,j,k)[IPV0_[aux1][aux2 + 1]];

                                    Xin = v(i,j,k)[IPV0_[aux1][aux2 + 1]];
                                }

                                e /= Foam::mag(e);

                                xgi(i,j,k) += Xin
                                        + (-C - (n(i,j,k) & Xin))
                                        * e * (1 / (n(i,j,k) & e));

                                numberOfPoints++;
                            }
                        }
                    }

                    xgi(i,j,k) /= numberOfPoints;

                    if (xgi(i,j,k)[1] != xgi(i,j,k)[1])
                    {
                        FatalErrorInFunction
                            <<  "LSFIR: center of vertex computation failed."
                            << endl << abort(FatalError);
                    }
                }

            }

        }

        xgi[0].correctBoundaryConditions();
        n[0].correctBoundaryConditions();

        // Minimize the square problem for all cells and store the solution

        forAllCells(n, i, j, k)
        {

            if
            (
                alpha_(i,j,k) > vof::threshold
             && alpha_(i,j,k) < 1 - vof::threshold
            )
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

                int nNeighs = 0;

                for (int aux1 = 0; aux1 < 3; aux1++)
                {
                    for (int aux2 = 0; aux2 < 3; aux2++)
                    {
                        for (int aux3 = 0; aux3 < 3; aux3++)
                        {
                            bool interiorNode =
                                (
                                    (i+aux1-1) >= n.I().left()
                                 && (i+aux1-1) <  n.I().right()
                                 && (j+aux2-1) >= n.I().bottom()
                                 && (j+aux2-1) <  n.I().top()
                                 && (k+aux3-1) >= n.I().aft()
                                 && (k+aux3-1) <  n.I().fore()
                                );

                            if
                            (
                                alpha_(i+aux1-1,j+aux2-1,k+aux3-1)
                              > vof::threshold
                             && alpha_(i+aux1-1,j+aux2-1,k+aux3-1)
                              < (1.0 - vof::threshold)
                             && (
                                    (aux1 != 1)
                                 || (aux2 != 1)
                                 || (aux3 != 1)
                                )
                             && (
                                    interiorNode
                                 || boundaryTypeLSFIR_[aux1][aux2][aux3]
                                )
                            )
                            {
                                scalar angle =
                                    (
                                        n(i+aux1-1,j+aux2-1,k+aux3-1)
                                      & n(i,j,k)
                                    )
                                  / (
                                        Foam::mag
                                        (
                                            n(i+aux1-1,j+aux2-1,k+aux3-1)
                                        )
                                      * Foam::mag(n(i,j,k))
                                    );

                                if (angle > validAngleTol)
                                {
                                    vector dist =
                                        xgi(i+aux1-1,j+aux2-1,k+aux3-1)
                                      - xgi(i,j,k);

                                    wi =
                                        1.0
                                      / Foam::max
                                        (
                                            1e-20,
                                            Foam::pow(Foam::mag(dist),2.5)
                                        );

                                    A[0][0] += wi*Foam::sqr(dist[d2]);
                                    A[0][1] += wi*(dist[d2]) * (dist[d3]);
                                    A[1][1] += wi*Foam::sqr(dist[d3]);

                                    b[0] += - wi*(dist[d2]) * (dist[d1]);
                                    b[1] += - wi*(dist[d3]) * (dist[d1]);

                                    nNeighs++;
                                }
                            }
                        }
                    }
                }

                A[1][0] = A[0][1];

                vector nc;

                if (nNeighs > 0)
                {
                    if ((A[0][0] * A[1][1] - A[0][1] * A[1][0]) == 0)
                    {
                        nc = n(i,j,k);
                    }
                    else
                    {
                        if
                        (
                            (A[0][0] == 0.0 && A[1][1] != 0)
                         || Foam::mag(A[0][0]/stabilise(A[1][1], 1e-40))
                          < 1e-12
                        )
                        {
                            nc[d2] = 0.0;
                            nc[d3] = b[1] / stabilise(A[1][1], 1e-40);
                        }
                        else if
                        (
                            (A[0][0] != 0.0 && A[1][1] == 0)
                         || Foam::mag(A[1][1]/stabilise(A[0][0], 1e-40))
                          < 1e-12
                        )
                        {
                            nc[d3] = 0.0;
                            nc[d2] = b[0] / stabilise(A[0][0], 1e-40);
                        }
                        else
                        {
                            nc[d2] =
                                (A[1][1] * b[0] - A[0][1] * b[1])
                              / (A[0][0] * A[1][1] - A[0][1] * A[1][0]);

                            nc[d3] =
                                (- A[1][0] * b[0] + A[0][0] * b[1])
                              / (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
                        }

                        if (n(i,j,k)[d1] < 0)
                        {
                            nc[d1] = - 1.0;
                            nc[d2] = - nc[d2];
                            nc[d3] = - nc[d3];
                        }
                        else
                        {
                            nc[d1] = 1.0;
                        }
                    }

                    nc /= Foam::mag(nc);
                    scalar angle =
                        (nc & n(i,j,k))
                      / (Foam::mag(nc) * Foam::mag(n(i,j,k)));

                    if (angle > maxAngleTol)
                    {
                        nNew(i,j,k) = nc;

                        if (angle < angleTol)
                            updatedNormals = 1;
                    }
                    else
                    {
                        nNew(i,j,k) = n(i,j,k);
                    }
                }
                else
                {
                    nNew(i,j,k) = n(i,j,k);
                }

                if
                (
                    nNew(i,j,k)[1] != nNew(i,j,k)[1]
                 || nNew(i,j,k)[2] != nNew(i,j,k)[2]
                 || nNew(i,j,k)[0] != nNew(i,j,k)[0]
                )
                {
                    FatalErrorInFunction
                        <<  "LSFIR: normal computation failed."
                        <<  "n new: " << nNew(i,j,k) << endl
                        << endl << abort(FatalError);
                }
            }
        }

        // Update the normals

        forAllCells(n, i, j, k)
        {
            if
            (
                alpha_(i,j,k) > vof::threshold
             && alpha_(i,j,k) < 1.0 - vof::threshold
            )
            {
                n(i,j,k) = nNew(i,j,k);
            }
        }
    }

    n[0].correctBoundaryConditions();
}

}

}

}
