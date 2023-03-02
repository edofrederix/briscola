#include "arcEdge.H"
#include "geometry.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(arcEdge, 0);
addToRunTimeSelectionTable(edge, arcEdge, dictionary);

arcEdge::arcEdge
(
    const face& f,
    const label num,
    const labelBlock& vertexNums,
    const labelVector& N
)
:
    edge(f, num, vertexNums, N),
    center_(),
    angle_(),
    axis_(),
    radius_(),
    start_()
{
    const vector p(dict().lookup("point"));

    const vector v0(this->v0());
    const vector v1(this->v1());

    const scalar pi = constant::mathematical::pi;

    // Relative to v0

    const vector a(p-v0);
    const vector b(v1-v0);

    const scalar aa(a&a);
    const scalar bb(b&b);
    const scalar ab(a&b);

    const scalar fact(0.5*(bb-ab)/(aa*bb-ab*ab));

    // Center of the circle (global)

    const vector c((0.5*a + fact*((a^b)^a)) + v0);

    // Relative to c

    const vector r0(v0-c);
    const vector rp(p-c);
    const vector r1(v1-c);

    scalar angle
    (
        Foam::acos
        (
            Foam::max
            (
                -1,
                Foam::min
                (
                    (r1 & r0)
                  / (mag(r1)*mag(r0)),
                    1
                )
            )
        )
    );

    if (((r0 ^ rp) & (r0 ^ r1)) < 0.0)
    {
        angle = 2.0*pi - angle;
    }

    vector axis;

    if (angle <= pi)
    {
        axis = r0^r1;
    }
    else
    {
        axis = r1^r0;
    }

    center_ = c;
    angle_ = angle;
    axis_ = axis/mag(axis);
    radius_ = mag(r0);
    start_ = r0/radius_;
}

arcEdge::arcEdge(const arcEdge& e)
:
    edge(e),
    center_(e.center_),
    angle_(e.angle_),
    axis_(e.axis_),
    radius_(e.radius_),
    start_(e.start_)
{}

tmp<vectorField> arcEdge::operator()(const scalarField& f) const
{
    const scalar pi = constant::mathematical::pi;

    tmp<vectorField> tp(new vectorField(f.size(), Zero));

    vectorField& p = tp.ref();

    // Draw the arc in the (x-y) plane, starting at y = 0 and rotating about the
    // z-axis with the computed angle

    p.replace(0, radius_*cos(angle_*f));
    p.replace(1, radius_*sin(angle_*f));

    // Compute the rotation tensor from z-axis to arc axis

    const vector z(0,0,1);
    const vector v(z ^ axis_);

    tensor R1;

    if (mag(v) < 1e-10)
    {
        if ((z & axis_) < 0)
        {
            // Rotate 180 degrees around the y-axis

            R1 = tensor(-1,0,0,0,1,0,0,0,-1);
        }
        else
        {
            // The z-axis and arc axis are the same

            R1 = I;
        }
    }
    else
    {
        // Rotation tensor from z to axis

        tensor W1
        (
            0, -v.z(), v.y(),
            v.z(), 0, -v.x(),
            -v.y(), v.x(), 0
        );

        R1 = I + W1 + (1.0-(z & axis_))/magSqr(v)*(W1 & W1);
    }

    // Computate the rotation tensor from rotated x-axis to the start of the
    // arc, about the arc axis. Note that this assumes that dot(axis,start) = 0

    const vector x(R1 & vector(1,0,0));

    tensor R2;

    scalar theta
    (
        Foam::acos
        (
            Foam::max
            (
                -1,
                Foam::min
                (
                    (x & start_)
                  / (mag(x)*mag(start_)),
                    1
                )
            )
        )
    );

    if ((axis_ & (x ^ start_)) < 0.0)
    {
        theta = 2.0*pi - theta;
    }

    const tensor W2
    (
        0, -axis_.z(), axis_.y(),
        axis_.z(), 0, -axis_.x(),
        -axis_.y(), axis_.x(), 0
    );

    R2 = I + Foam::sin(theta)*W2 + 2.0*sqr(Foam::sin(theta/2.0))*(W2 & W2);

    p = (R2 & (R1 & p)) + center_;

    return tp;
}

}

}
