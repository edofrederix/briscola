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

using Foam::constant::mathematical::pi;

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
    LSGIR(fvMsh, dict, alpha),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 2)),
    tolerance_
    (
        dict.lookupOrDefault<scalar>
        (
            "tolerance",
            2.0*pi/360.0
        )
    )
{}

LSFIR::LSFIR(const LSFIR& s)
:
    LSGIR(s),
    maxIter_(s.maxIter_),
    tolerance_(s.tolerance_)
{}

LSFIR::~LSFIR()
{}

void LSFIR::correct()
{
    // Reconstruct the interface normal using the version of LSGIR presented in
    // Lopez & Hernández (2022)

    LSGIR::correct();

    colocatedVectorField& n = *this;

    const faceLabel I(n.I());

    const colocatedVertexVectorField& v =
        fvMsh_.template metrics<colocated>().vertexCenters();

    colocatedVectorField centroids("centroids", fvMsh_);

    // Temporary normal field constructed from the LSGIR normal

    colocatedVectorField nt(n);

    const LVE lve(fvMsh_.msh().rectilinear() == unitXYZ);

    label iter;
    for (iter = 0; iter < maxIter_; iter++)
    {
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

        const scalar maxChangeInAngle = 0.166*pi/scalar(iter+1);
        bool converged = true;

        forAllCells(nt, i, j, k)
        {
            const labelVector ijk(i,j,k);

            if
            (
                alpha_(ijk) > vof::threshold
             && alpha_(ijk) < 1 - vof::threshold
            )
            {
                // Find dominant normal index

                label d = 0;
                if (Foam::mag(n(ijk)[1]) > Foam::mag(n(ijk)[d])) d = 1;
                if (Foam::mag(n(ijk)[2]) > Foam::mag(n(ijk)[d])) d = 2;

                // Set other two indices (rotate)

                const label d1 = (d + 1)%3;
                const label d2 = (d + 2)%3;

                tensor2D A(tensor2D::zero);
                vector2D b(vector2D::zero);

                label count = 0;

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
                        // Only use neighbor if the angle between the normals is
                        // within 45 degrees

                        if (Foam::acos(n(ijk+o) & n(ijk)) < pi/4.0)
                        {
                            const vector dist =
                                centroids(ijk+o) - centroids(ijk);

                            const scalar weight =
                                1.0/Foam::pow(Foam::max(mag(dist), 1e-20), 2.5);

                            A.xx() += weight*Foam::sqr(dist[d1]);
                            A.xy() += weight*dist[d1]*dist[d2];
                            A.yy() += weight*Foam::sqr(dist[d2]);

                            b.x() -= weight*dist[d1]*dist[d];
                            b.y() -= weight*dist[d2]*dist[d];

                            count++;
                        }
                    }
                }

                A.yx() = A.xy();

                if (count > 0 && det(A) != 0.0)
                {
                    const scalar r = Foam::mag(A.xx()/stabilise(A.yy(), 1e-40));

                    vector nNew(Zero);

                    if (r < 1e-12)
                    {
                        nNew[d2] = b.y()/stabilise(A.yy(), 1e-40);
                    }
                    else if (1.0/r < 1e-12)
                    {
                        nNew[d1] = b.x()/stabilise(A.xx(), 1e-40);
                    }
                    else
                    {
                        const vector2D x(inv(A) & b);

                        nNew[d1] = x.x();
                        nNew[d2] = x.y();
                    }

                    nNew[d] = 1.0;
                    nNew /= Foam::mag(nNew);

                    if (n(ijk)[d] < 0)
                        nNew = -nNew;

                    const scalar angle = Foam::acos(nNew & n(ijk));

                    // Update only if the change in angle is smaller than a
                    // threshold

                    if (angle < maxChangeInAngle)
                    {
                        nt(ijk) = nNew;

                        // If the change in angle is larger than the angle
                        // tolerance, keep iterating

                        if (angle > tolerance_)
                            converged = false;
                    }
                }
            }
        }

        // Copy back the improved normals

        n = nt;

        reduce(converged, andOp<Switch>());

        if (converged)
            break;

        iter++;
    }

    n.correctBoundaryConditions();
}

}

}

}
