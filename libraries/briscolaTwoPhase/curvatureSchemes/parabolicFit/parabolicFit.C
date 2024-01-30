#include "parabolicFit.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"
#include "rectilinearMesh.H"
#include "LVE.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(parabolicFit, 0);
addToRunTimeSelectionTable(curvatureScheme, parabolicFit, dictionary);

void parabolicFit::createBoundaryTypes()
{
    const faceLabel& faceType = fvMsh_.msh().faceBoundaryType();

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
                        i == 0
                      ? aux2 && faceType.left()
                      > domainBoundary::typeNumber
                      : aux2 && faceType.right()
                      > domainBoundary::typeNumber;
                }

                if (j != 1)
                {
                    aux2 =
                        j == 0
                      ? aux2 && faceType.bottom()
                      > domainBoundary::typeNumber
                      : aux2 && faceType.top()
                      > domainBoundary::typeNumber;
                }

                if (k != 1)
                {
                    aux2 =
                        k == 0
                      ? aux2 && faceType.aft()
                      > domainBoundary::typeNumber
                      : aux2 && faceType.fore()
                      > domainBoundary::typeNumber;
                }

                boundaryType_[i][j][k] = aux2;
            }
        }
    }
}

void parabolicFit::createHexagonDescription()
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

parabolicFit::parabolicFit
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
:
    curvatureScheme(fvMsh, dict, normal, alpha),
    lve_(fvMsh.msh()[0].rectilinear() == unitXYZ)
{
    createBoundaryTypes();
    createHexagonDescription();
}

parabolicFit::parabolicFit(const parabolicFit& s)
:
    curvatureScheme(s),
    lve_(fvMsh_.msh()[0].rectilinear() == unitXYZ)
{
    createBoundaryTypes();
    createHexagonDescription();
}


parabolicFit::~parabolicFit()
{}

void parabolicFit::correct()
{

    colocatedScalarField& kappa = *this;

    const colocatedVertexVectorField& v =
        fvMsh_.template metrics<colocated>().vertexCenters();
    const colocatedScalarField& cv =
        fvMsh_.template metrics<colocated>().cellVolumes();

    tmp<colocatedVectorField> txgi
    (
        new colocatedVectorField
        (
            "surface_centers",
            fvMsh_
        )
    );

    colocatedVectorField& xgi = txgi.ref();

    forAllCells(kappa, i, j, k)
    {
        xgi(i,j,k) = Zero;

        if
        (
            alpha_(i,j,k) > vof::threshold
         && alpha_(i,j,k) < (1.0 - vof::threshold)
        )
        {

            // Prevent computing interface for zero normal

            //onst scalar S = Foam::mag(normal_(i,j,k));

            //if (S < 1e-12)
            //{
            //    n(i,j,k)[0] = 1;
            //}

            // Compute which vertex are inside the liquid for the truncated
            // cell

            scalar C =
                lve_(alpha_(i,j,k),v(i,j,k),cv(i,j,k),normal_(i,j,k));

            double originalCs[8];

            for (int aux1 = 0; aux1 < 8; aux1++)
            {
                originalCs[aux1] = (normal_(i,j,k) & v(i,j,k)[aux1]);
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


            int NIPV1[16] = {0};
            int IPV1[16][14];

            int insertedVertex[12][6];
            int totalInsertedVertex = 0;
            int Nf = originalNf_;
            int Ip = 8;

            for (int jt = 0; jt < originalNf_; jt++)
            {
                for (int it = 1; it <= NIPV0_[jt]; it++)
                    if (IA[IPV0_[jt][it]] == 1)
                        NIPV1[jt]++;

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
                                    insertedVertex[totalInsertedVertex][3]
                                        = jt;
                                    insertedVertex[totalInsertedVertex][4]
                                        = index;
                                    insertedVertex[totalInsertedVertex][2]
                                        = ip1;
                                    insertedVertex[totalInsertedVertex][1]
                                        = ip2;
                                }
                                else
                                {
                                    insertedVertex[totalInsertedVertex][1]
                                        = ip1;
                                    insertedVertex[totalInsertedVertex][2]
                                        = ip2;
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
                        << "Curvature (parabolicFit): Arrangement of the new face failed."
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
                    x0[insertedVertex[it][0] - 8] = Xin
                            + (-C - (normal_(i,j,k) & Xin))
                                * e * (1 / (normal_(i,j,k) & e));
                }

                scalar auxMax;
                int maxNormalIndex;

                auxMax = Foam::mag(normal_(i,j,k)[0]);
                maxNormalIndex = 0;

                for (int it = 1; it < 3; it++)
                {
                    if (Foam::mag(normal_(i,j,k)[it]) > auxMax)
                    {
                        auxMax = Foam::mag(normal_(i,j,k)[it]);
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

                    Area /= normal_(i,j,k)[maxNormalIndex];

                    TotalArea += Area;
                    xgi(i,j,k) += Area * Centroid;
                }

                xgi(i,j,k) /= TotalArea;

                if (xgi(i,j,k)[1] != xgi(i,j,k)[1])
                {
                    FatalErrorInFunction
                        << "Curvature (parabolicFit): geometric center computation failed."
                        << normal_(i,j,k) << endl
                        << C << endl
                        << xgi(i,j,k) << endl
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

                            xgi(i,j,k) +=
                                Xin
                              + (-C - (normal_(i,j,k) & Xin))
                              * e * (1 / (normal_(i,j,k) & e));

                            numberOfPoints++;
                        }
                    }
                }

                xgi(i,j,k) /= numberOfPoints;

                if (xgi(i,j,k)[1] != xgi(i,j,k)[1])
                {
                    FatalErrorInFunction
                        << "Curvature (parabolicFit): center of vertex computation failed."
                        << endl << abort(FatalError);
                }
            }
        }
    }

    txgi.ref()[0].correctBoundaryConditions();

    forAllCells(kappa, i, j, k)
    {

        if
        (
            alpha_(i,j,k) >       vof::threshold
         && alpha_(i,j,k) < 1.0 - vof::threshold
        )
        {
            if (fvMsh_[0].n() > 1)
            {
                double A[6][6] = {0};
                double L[6][6] = {0};
                double y[6];
                double coefs[6];
                double b[6] = {0};
                vector auxBase(1,0,0);

                if (Foam::mag(auxBase ^ normal_(i,j,k)) < 1e-2)
                    auxBase = vector(0,1,0);

                vector v1 = normal_(i,j,k) ^ auxBase;
                v1 /= Foam::mag(v1);
                vector v2 = normal_(i,j,k) ^ v1;
                v2 /= Foam::mag(v2);

                const tensor T
                (
                    v1,
                    v2,
                    normal_(i,j,k) / Foam::mag(normal_(i,j,k))
                );

                int count = 0;

                for (int aux1 = -1; aux1 <= 1; aux1++)
                {
                    for (int aux2 = -1; aux2 <= 1; aux2++)
                    {
                        for (int aux3 = -1; aux3 <= 1; aux3++)
                        {
                            bool interiorNode =
                                (
                                    (i+aux1-1) >= kappa.I().left()
                                 && (i+aux1-1) <  kappa.I().right()
                                 && (j+aux2-1) >= kappa.I().bottom()
                                 && (j+aux2-1) <  kappa.I().top()
                                 && (k+aux3-1) >= kappa.I().aft()
                                 && (k+aux3-1) <  kappa.I().fore()
                                );

                            if
                            (
                                alpha_(i+aux1,j+aux2,k+aux3) > vof::threshold
                             && alpha_(i+aux1,j+aux2,k+aux3)
                              < 1.0 - vof::threshold
                             && (
                                    interiorNode
                                 || boundaryType_[aux1 + 1][aux2 + 1][aux3 + 1]
                                )
                            )
                            {
                                vector xNew =
                                    T
                                  & (xgi(i+aux1,j+aux2,k+aux3) - xgi(i,j,k));

                                A[0][0] += Foam::pow(xNew.x(),4);
                                A[0][1] +=
                                    Foam::pow(xNew.x(),2)
                                  * Foam::pow(xNew.y(),2);
                                A[0][2] +=
                                    Foam::pow(xNew.x(),3)
                                  * Foam::pow(xNew.y(),1);
                                A[0][3] += Foam::pow(xNew.x(),3);
                                A[0][4] +=
                                    Foam::pow(xNew.x(),2)
                                  * Foam::pow(xNew.y(),1);
                                A[0][5] += Foam::pow(xNew.x(),2);

                                A[1][1] += Foam::pow(xNew.y(),4);
                                A[1][2] +=
                                    Foam::pow(xNew.x(),1)
                                  * Foam::pow(xNew.y(),3);
                                A[1][3] +=
                                    Foam::pow(xNew.x(),1)
                                  * Foam::pow(xNew.y(),2);
                                A[1][4] += Foam::pow(xNew.y(),3);
                                A[1][5] += Foam::pow(xNew.y(),2);

                                A[2][2] +=
                                    Foam::pow(xNew.x(),2)
                                  * Foam::pow(xNew.y(),2);
                                A[2][3] +=
                                    Foam::pow(xNew.x(),2)
                                  * Foam::pow(xNew.y(),1);
                                A[2][4] +=
                                    Foam::pow(xNew.x(),1)
                                  * Foam::pow(xNew.y(),2);
                                A[2][5] +=
                                    Foam::pow(xNew.x(),1)
                                  * Foam::pow(xNew.y(),1);

                                A[3][3] += Foam::pow(xNew.x(),2);
                                A[3][4] +=
                                    Foam::pow(xNew.x(),1)
                                  * Foam::pow(xNew.y(),1);
                                A[3][5] += Foam::pow(xNew.x(),1);

                                A[4][4] += Foam::pow(xNew.y(),2);
                                A[4][5] += Foam::pow(xNew.y(),1);

                                A[5][5] += 1;

                                b[0] += xNew.z() * Foam::pow(xNew.x(),2);
                                b[1] += xNew.z() * Foam::pow(xNew.y(),2);
                                b[2] +=
                                    xNew.z() * Foam::pow(xNew.x(),1)
                                  * Foam::pow(xNew.y(),1);
                                b[3] += xNew.z() * Foam::pow(xNew.x(),1);
                                b[4] += xNew.z() * Foam::pow(xNew.y(),1);
                                b[5] += xNew.z();

                                count++;
                            }
                        }
                    }
                }

                if (count > 5)
                {

                    for (int it = 1; it < 6; it++)
                        for(int jt = 0; jt < it; jt++)
                            A[it][jt] = A[jt][it];

                    for (int it = 0; it < 6; it++)
                    {
                        for (int jt = 0; jt < 6; jt++)
                        {
                            if (jt >= it)
                            {
                                L[jt][it] = A[jt][it];
                                for (int kt = 0; kt < it; kt++)
                                {
                                    L[jt][it] = L[jt][it] - L[jt][kt] * L[kt][it];
                                }
                            }
                        }
                        for (int jt = 0; jt < 6; jt++)
                        {
                            if (jt > it)
                            {
                                L[it][jt] = A[it][jt] / L[it][it];
                                for (int kt = 0; kt < it; kt++)
                                {
                                    L[it][jt] = L[it][jt] - ((L[it][kt] * L[kt][jt]) / L[it][it]);
                                }
                            }
                        }
                    }

                    for(int it = 0; it < 6; it++)
                    {
                        double sum = 0;

                        for(int jt = 0; jt < it; jt++)
                            sum = sum + L[it][jt]*y[jt];

                        y[it] = (b[it] - sum)/L[it][it];
                    }

                    for(int it = 5; it >= 0; it--)
                    {
                        double sum = 0;

                        for(int jt = it + 1; jt < 6; jt++)
                            sum = sum + L[it][jt]*coefs[jt];

                        coefs[it] = y[it] - sum;
                    }

                    kappa(i,j,k) = 2 *
                            (
                                coefs[0] * (1 + Foam::sqr(coefs[4]))
                                + coefs[1] * (1 + Foam::sqr(coefs[3]))
                                - coefs[2] * coefs[3] * coefs[4]
                            )
                          / (
                                Foam::pow
                                (
                                    1 + Foam::sqr(coefs[3])
                                  + Foam::sqr(coefs[4]),
                                    1.5
                                )
                            );

                    scalar maxKappa = 1/Foam::pow(cv(i,j,k),1.0/3.0);

                    if(Foam::mag(kappa(i,j,k)) > maxKappa)
                    {
                        kappa(i,j,k) = maxKappa * Foam::sign(kappa(i,j,k));
                    }

                    if (kappa(i,j,k) != kappa(i,j,k))
                    {
                        FatalErrorInFunction
                            << "Curvature (parabolicFit): kappa computation failed."
                            << endl << abort(FatalError);
                    }
                }
                else
                {
                    kappa(i,j,k) = 0;
                }
            }
            else
            {
                double A[3][3] = {0};
                double L[3][3] = {0};
                double y[3];
                double coefs[3];
                double b[3] = {0};

                int count = 0;

                const tensor T
                (
                    vector(normal_(i,j,k)[1], -normal_(i,j,k)[0], 0),
                    vector(normal_(i,j,k)[0], normal_(i,j,k)[1], 0),
                    vector(0, 0, 1)
                );

                for (int aux1 = -1; aux1 <= 1; aux1++)
                {
                    for (int aux2 = -1; aux2 <= 1; aux2++)
                    {
                        for (int aux3 = -1; aux3 <= 1; aux3++)
                        {
                            bool interiorNode =
                                    (
                                        (i+aux1) >= kappa.I().left()
                                     && (i+aux1) <  kappa.I().right()
                                     && (j+aux2) >= kappa.I().bottom()
                                     && (j+aux2) <  kappa.I().top()
                                     && (k+aux3) >= kappa.I().aft()
                                     && (k+aux3) <  kappa.I().fore()
                                    );

                            if
                            (
                                alpha_(i+aux1,j+aux2,k+aux3) > vof::threshold
                             && alpha_(i+aux1,j+aux2,k+aux3)
                              < (1.0 - vof::threshold)
                             && (
                                    interiorNode
                                 || boundaryType_[aux1 + 1][aux2 + 1][aux3 + 1]
                                )
                            )
                            {
                                vector xNew =
                                    T
                                  & (xgi(i+aux1,j+aux2,k+aux3) - xgi(i,j,k));

                                A[0][0] += Foam::pow(xNew.x(),4);
                                A[0][1] += Foam::pow(xNew.x(),3);
                                A[0][2] += Foam::pow(xNew.x(),2);

                                A[1][1] += Foam::pow(xNew.x(),2);
                                A[1][2] += xNew.x();

                                A[2][2] += 1;

                                b[0] += xNew.y() * Foam::pow(xNew.x(),2);
                                b[1] += xNew.y() * xNew.x();
                                b[2] += xNew.y();

                                count++;

                            }
                        }
                    }
                }

                for (int it = 1; it < 3; it++)
                    for(int jt = 0; jt < it; jt++)
                        A[it][jt] = A[jt][it];

                if (count > 2)
                {
                    for (int it = 0; it < 3; it++)
                    {
                        for (int jt = 0; jt < 3; jt++)
                        {
                            if (jt >= it)
                            {
                                L[jt][it] = A[jt][it];
                                for (int kt = 0; kt < it; kt++)
                                {
                                    L[jt][it] = L[jt][it] - L[jt][kt] * L[kt][it];
                                }
                            }
                        }
                        for (int jt = 0; jt < 3; jt++)
                        {
                            if (jt > it)
                            {
                                L[it][jt] = A[it][jt] / L[it][it];
                                for (int kt = 0; kt < it; kt++)
                                {
                                    L[it][jt] = L[it][jt] - ((L[it][kt] * L[kt][jt]) / L[it][it]);
                                }
                            }
                        }
                    }

                    for(int it = 0; it < 3; it++)
                    {
                        double sum = 0;

                        for(int jt = 0; jt < it; jt++)
                            sum = sum + L[it][jt]*y[jt];

                        y[it] = (b[it] - sum)/L[it][it];
                    }

                    for(int it = 2; it >= 0; it--)
                    {
                        double sum = 0;

                        for(int jt = it + 1; jt < 3; jt++)
                            sum = sum + L[it][jt]*coefs[jt];

                        coefs[it] = y[it] - sum;
                    }

                    kappa(i,j,k) = 2 *
                            (
                                coefs[0]
                            )
                          / (
                                Foam::pow(1 + Foam::sqr(coefs[1]), 1.5)
                            );

                    scalar maxKappa = 1/Foam::pow(cv(i,j,k)/0.1,1.0/2.0);

                    if(Foam::mag(kappa(i,j,k)) > maxKappa)
                    {
                        kappa(i,j,k) = maxKappa * Foam::sign(kappa(i,j,k));
                    }


                }
                else
                {
                    kappa(i,j,k) = 0;
                }
            }
        }
        else
        {
            kappa(i,j,k) = 0;
        }
    }

    kappa.correctBoundaryConditions();
}

}

}

}
