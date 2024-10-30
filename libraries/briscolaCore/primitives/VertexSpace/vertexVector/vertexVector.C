#include "vertexVector.H"

namespace Foam
{

template<>
const char* const vertexVector::csType::typeName = "vertexVector";

template<>
const char* const vertexVector::csType::componentNames[] =
    {"v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7"};

template<>
const labelVector vertexVector::csType::componentOffsets[] =
    {
        briscola::vertexOffset0,
        briscola::vertexOffset1,
        briscola::vertexOffset2,
        briscola::vertexOffset3,
        briscola::vertexOffset4,
        briscola::vertexOffset5,
        briscola::vertexOffset6,
        briscola::vertexOffset7,
    };

template<>
const vertexVector vertexVector::csType::zero
(
    vertexVector::uniform(vector::uniform(0))
);

template<>
const vertexVector vertexVector::csType::one
(
    vertexVector::uniform(vector::uniform(1))
);

template<>
const vertexVector vertexVector::csType::max
(
    vertexVector::uniform(vector::uniform(vGreat))
);

template<>
const vertexVector vertexVector::csType::min
(
    vertexVector::uniform(vector::uniform(-vGreat))
);

template<>
const vertexVector vertexVector::csType::rootMax
(
    vertexVector::uniform(vector::uniform(rootVGreat))
);

template<>
const vertexVector vertexVector::csType::rootMin
(
    vertexVector::uniform(vector::uniform(-rootVGreat))
);

vector interpolationWeights
(
    const vector& point,
    const vertexVector& vertices,
    const bool insideOnly,
    const bool fatal
)
{
    const scalar tol = 1e-14;

    if (insideOnly)
    {
        // Bounding box check for early rejection

        scalar xMin(vertices[0].x());
        scalar xMax(vertices[0].x());

        scalar yMin(vertices[0].y());
        scalar yMax(vertices[0].y());

        scalar zMin(vertices[0].z());
        scalar zMax(vertices[0].z());

        for (int ii = 1; ii < 8; ii++)
        {
            xMin = Foam::min(xMin,vertices[ii].x());
            xMax = Foam::max(xMax,vertices[ii].x());

            yMin = Foam::min(yMin,vertices[ii].y());
            yMax = Foam::max(yMax,vertices[ii].y());

            zMin = Foam::min(zMin,vertices[ii].z());
            zMax = Foam::max(zMax,vertices[ii].z());
        }

        if
        (
            point.x() < xMin-tol || point.x() > xMax+tol
         || point.y() < yMin-tol || point.y() > yMax+tol
         || point.z() < zMin-tol || point.z() > zMax+tol
        )
        {
            return -vector::one;
        }
    }

    // Try to find the parameters (u,v,w) in the parametrized function that
    // tri-linearly interpolates the hexahedron. This function properly accounts
    // for faces being doubly ruled surfaces. If (u,v,w) is in the unit cube,
    // the point is in the cell. Find (u,v,w) with Newton's method.

    const vector P(point - vertices[0]);
    const vertexVector D(vertices - vertices[0]);

    const vector F1(D.rba());
    const vector F2(D.lta());
    const vector F3(D.lbf());

    const vector F4(D.rta()-F1-F2);
    const vector F5(D.rbf()-F1-F3);
    const vector F6(D.ltf()-F2-F3);

    const vector F7(D.rtf()-F1-F2-F3-F4-F5-F6);

    // Initial guess (exact for rectangular cells)

    vector u
    (
        (P & D.rba())/Foam::magSqr(D.rba()),
        (P & D.lta())/Foam::magSqr(D.lta()),
        (P & D.lbf())/Foam::magSqr(D.lbf())
    );

    const label maxIter = 100;

    vector du(vector::one);
    label iter = 0;

    while (mag(du) > tol && iter < maxIter)
    {
        const tensor dfdu
        (
            F1 + u.y()*F4 + u.z()*F5 + u.y()*u.z()*F7,
            F2 + u.x()*F4 + u.z()*F6 + u.x()*u.z()*F7,
            F3 + u.x()*F5 + u.y()*F6 + u.x()*u.y()*F7
        );

        const vector f
        (
            u.x()*F1
          + u.y()*F2
          + u.z()*F3
          + u.x()*u.y()*F4
          + u.x()*u.z()*F5
          + u.y()*u.z()*F6
          + u.x()*u.y()*u.z()*F7
          - P
        );

        du = inv(dfdu.T()) & f;
        u -= du;

        iter++;
    }

    if (iter == maxIter)
    {
        if (fatal)
            FatalErrorInFunction
                << "Could not determine the interpolation weights of point "
                << point << " in the hexahedron " << vertices << endl
                << abort(FatalError);

        return -vector::one;
    }

    // Round up to tol. This seems to be important to have the point belong to
    // the correct cell when the point is very close to a face.

    u = vector
    (
        round(u.x()/tol)*tol,
        round(u.y()/tol)*tol,
        round(u.z()/tol)*tol
    );

    if
    (
        // At faces points belong to the upper cell

        insideOnly
     && (
            u.x() < 0 || u.x() >= 1
         || u.y() < 0 || u.y() >= 1
         || u.z() < 0 || u.z() >= 1
        )
    )
    {
        return -vector::one;
    }
    else
    {
        return u;
    }
}

}
