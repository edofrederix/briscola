#include "LSFIR.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"
#include "LVE.H"
#include "SortableList.H"
#include "constants.H"
#include "Time.H"
#include "tensor2D.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(LSFIR, 0);
addToRunTimeSelectionTable(normalScheme, LSFIR, dictionary);

vector LSFIR::centroid
(
    const labelVector& ijk,
    const List<Switch>& IA,
    const scalar C
) const
{
    const vertexVector v =
        fvMsh_.template metrics<colocated>().vertexCenters()(ijk);

    const vector n = this->operator()(ijk);

    // Collect intersection points in unsorted list

    DynamicList<vector> uPoints;

    for (label fi = 0; fi < 12; fi++)
    for (label vi = 0; vi < 3; vi++)
    {
        const label v1 = faceTriangleDecompCC[fi][vi];
        const label v2 = faceTriangleDecompCC[fi][(vi+1)%3];

        if (IA[v1] != IA[v2])
        {
            const vector x1 = v[v1];
            const vector x2 = v[v2];

            const vector d(x2 - x1);

            if (Foam::mag(n & d) > VSMALL)
                uPoints.append(trimPrecision(x1 - (C + (n & x1))*d/(n & d)));
        }
    }

    if (uPoints.empty())
        return Zero;

    // Find dominant normal index

    label d = 0;
    if (Foam::mag(n[1]) > Foam::mag(n[d])) d = 1;
    if (Foam::mag(n[2]) > Foam::mag(n[d])) d = 2;

    // Set other two indices (rotate)

    const label d1 = (d + 1)%3;
    const label d2 = (d + 2)%3;

    // Compute 2D centroid as the average of all points

    vector2D average(vector2D::zero);

    forAll(uPoints, i)
    {
        average.x() += uPoints[i][d1];
        average.y() += uPoints[i][d2];
    }

    average /= scalar(uPoints.size());

    // Compute angles around the centroid

    scalarList angles(uPoints.size());
    forAll(uPoints, i)
    {
        const scalar ax = uPoints[i][d1] - average.x();
        const scalar ay = uPoints[i][d2] - average.y();

        angles[i] = Foam::atan2(ay,ax);
    }

    // Sort unique points according to their angles

    labelList order;
    uniqueOrder(angles, order);

    vectorList points(order.size());
    forAll(points, i)
        points[i] = uPoints[order[i]];

    // Compute centroid of polygon

    const vector p0 = points[0];
    scalar W = 0.0;
    vector centroid(Zero);

    for (label i = 1; i < points.size()-1; i++)
    {
        const vector a1 = points[i] - p0;
        const vector a2 = points[i+1] - p0;

        const scalar w = Foam::mag(a1 ^ a2)/2.0;
        const vector c = p0 + (a1 + a2)/3.0;

        W += w;
        centroid += w*c;
    }

    return
        Foam::mag(W) > VSMALL
      ? trimPrecision(centroid/W)
      : centroidFromVertices(ijk, IA, C);
}

vector LSFIR::centroidFromVertices
(
    const labelVector& ijk,
    const List<Switch>& IA,
    const scalar C
) const
{
    const vertexVector v =
        fvMsh_.template metrics<colocated>().vertexCenters()(ijk);

    const vector n = this->operator()(ijk);

    vector c(Zero);

    label i = 0;
    for (label fi = 0; fi < 12; fi++)
    for (label vi = 0; vi < 3;  vi++)
    {
        const label v1 = faceTriangleDecompCC[fi][vi];
        const label v2 = faceTriangleDecompCC[fi][(vi+1)%3];

        if (IA[v1] != IA[v2])
        {
            const vector x1 = v[v1];
            const vector x2 = v[v2];

            const vector d(x2 - x1);

            if (Foam::mag(n & d) > VSMALL)
            {
                c += x1 - (C + (n & x1))*d/(n & d);
                i++;
            }
        }
    }

    if (i == 0)
        return Zero;

    return trimPrecision(c/scalar(i));
}

LSFIR::LSFIR
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
:
    LSGIR(fvMsh, dict, alpha)
{}

LSFIR::LSFIR(const LSFIR& s)
:
    LSGIR(s)
{}

LSFIR::~LSFIR()
{}

void LSFIR::correct()
{
    LSGIR::correct();

    const scalar pi = Foam::constant::mathematical::pi;

    colocatedVectorField& n = *this;

    const faceLabel I(n.I());

    // Reconstruct the interface normal using the version of LSGIR presented in
    // Lopez & Hernández (2022)

    const colocatedVertexVectorField& v =
        fvMsh_.template metrics<colocated>().vertexCenters();

    colocatedVectorField centroids("centroids", fvMsh_);

    const scalar angleTol = Foam::cos(1e-3);
    const scalar validAngleTol = Foam::cos(pi/4);

    // Temporary normal field constructed from the LSGIR normal

    colocatedVectorField nt(n);

    const LVE lve(fvMsh_.msh().rectilinear() == unitXYZ);

    Switch update = true;
    label iter = 0;

    while (update && iter < 1)
    {
        const scalar maxAngleTol =
            Foam::cos(0.166*pi/scalar(iter+1));

        update = false;

        // Compute the centroids of the interface for all interfacial cells

        centroids = Zero;

        forAllCells(centroids, i, j, k)
        {
            const labelVector ijk(i,j,k);

            if
            (
                alpha_(ijk) > vof::threshold
             && alpha_(ijk) < 1- vof::threshold
            )
            {
                // Compute the C coefficient based on the LSGIR normal

                scalar Ci =
                    trimPrecision(lve(alpha_(ijk), v(ijk), n(ijk)));

                // Compute vertex C coefficients

                scalarList C(8);
                forAll(C, ii)
                    C[ii] = trimPrecision(n(ijk) & v(ijk)[ii]);

                // Check if vertex C values are larger than Ci

                List<Switch> IA(8);
                forAll(IA, ii)
                    IA[ii] = C[ii] < -Ci;

                // Compute the plane centroid

                centroids(ijk) = centroid(ijk, IA, Ci);
            }
        }

        centroids.correctBoundaryConditions();
        n.correctBoundaryConditions();

        // Minimize the square problem for all cells and store the solution

        forAllCells(nt, i, j, k)
        {
            const labelVector ijk(i,j,k);

            if
            (
                alpha_(ijk) > vof::threshold
             && alpha_(ijk) < 1 - vof::threshold
            )
            {
                // Get order of normal component magnitudes

                labelList order;
                sortedOrder(list(cmptMag(n(ijk))), order);

                // Collect 2x2 linear system

                label count = 0;
                tensor2D A(tensor2D::zero);
                vector2D b(vector2D::zero);

                labelVector o;
                for (o.x() = -1; o.x() < 2; o.x()++)
                for (o.y() = -1; o.y() < 2; o.y()++)
                for (o.z() = -1; o.z() < 2; o.z()++)
                if (o != zeroXYZ)
                {
                    const bool internal =
                        cmptMin(ijk + o - I.lower()) > -1
                     && cmptMax(ijk + o - I.upper()) <  0;

                    const bool pBoundary =
                        fvMsh_.msh().pBoundaryMask()(o + unitXYZ);

                    const bool interface =
                        alpha_(ijk + o) > vof::threshold
                     && alpha_(ijk + o) < 1 - vof::threshold;

                    if (interface && (internal || pBoundary))
                    {
                        const scalar angle = (n(ijk+o) & n(ijk));

                        Pout<< n(ijk) << " " << n(ijk+o) << endl;

                        if (angle > validAngleTol)
                        {
                            const vector dist =
                                centroids(ijk+o) - centroids(ijk);

                            const scalar weight =
                                1.0/Foam::pow(Foam::max(mag(dist), 1e-16), 2.5);

                            A.xx() += weight*Foam::sqr(dist[order[0]]);
                            A.yy() += weight*Foam::sqr(dist[order[1]]);

                            A.xy() += weight*dist[order[0]]*dist[order[1]];
                            A.yx() += weight*dist[order[0]]*dist[order[1]];

                            b.x() -= weight*dist[order[2]]*dist[order[0]];
                            b.y() -= weight*dist[order[2]]*dist[order[1]];

                            count++;
                        }
                    }
                }

                // Remove truncation errors

                for (int ii = 0; ii < 4; ii++)
                    A[ii] = trimPrecision(A[ii]);

                // Solve 2x2 linear system if a solution exists

                if (count > 0 && det(A) != 0.0)
                {
                    vector2D x(inv(A) & b);

                    vector nNew;

                    nNew[order[0]] = x.x();
                    nNew[order[1]] = x.y();
                    nNew[order[2]] = 1.0;

                    if (n(ijk)[order[0]] < 0)
                        nNew = -nNew;

                    nNew /= mag(nNew);

                    const scalar angle = (nNew & n(ijk));

                    if (angle > maxAngleTol)
                    {
                        nt(ijk) = nNew;

                        if (angle < angleTol)
                            update = true;
                    }
                }
            }
        }

        // Copy back the improved normals

        n = nt;

        reduce(update, orOp<Switch>());

        iter++;
    }

    n.correctBoundaryConditions();
}

}

}

}
