#include "truncatedPiped.H"

namespace Foam
{

namespace briscola
{

namespace fv
{


truncatedPiped::truncatedPiped
(
    const vector& v_lba,
    const vector& v_rba,
    const vector& v_lta,
    const vector& v_lbf,
    const vector& n,
    const scalar& C
)
:
    v_lba_(v_lba),
    v_rba_(v_rba),
    v_lta_(v_lta),
    v_lbf_(v_lbf),
    n_(n),
    C_(C)
{}


truncatedPiped::truncatedPiped
(
    const vertexVector& v,
    const vector& n,
    const scalar& C
)
:
    v_lba_(v.lba()),
    v_rba_(v.rba()),
    v_lta_(v.lta()),
    v_lbf_(v.lbf()),
    n_(n),
    C_(C)
{}

truncatedPiped::truncatedPiped(const truncatedPiped& hx)
:
    v_lba_(hx.v_lba_),
    v_rba_(hx.v_rba_),
    v_lta_(hx.v_lta_),
    v_lbf_(hx.v_lbf_),
    n_(hx.n_),
    C_(hx.C_)
{}

truncatedPiped::~truncatedPiped()
{}

scalar truncatedPiped::volume() const
{

    /*

    This function computes and returns the volume of the parallelepiped
    truncated by the plane x * n + C = 0.

    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    First the linear transform

    x_old = v_.lba + T * x_new

    projects the parallelepiped into the unit cube, then the algorithm in
    Scardovelli & Zaleski (2000) is used.

    */

    const tensor T
    (
        v_rba_ - v_lba_,
        v_lta_ - v_lba_,
        v_lbf_ - v_lba_
    );

    vector n = T & n_;
    scalar C = C_ + (n_ & v_lba_);

    for (int i = 0; i < 3; i++)
    {
       if (n[i] < 0)
        {
            C += n[i];
            n[i] = -n[i];
        }
    }

    scalar scalingFactor = Foam::mag(Foam::det(T));

    scalar s = cmptSum(cmptMag(n));
    n = cmptMag(n)/s;
    C = -C/s;

    scalarList ml(3);

    for (int i = 0; i < 3; i++)
        ml[i] = n[i];

    sort(ml);

    scalar& m1 = ml[0];
    scalar& m2 = ml[1];
    scalar& m3 = ml[2];

    scalar m12 = m1 + m2;
    scalar V, mm;
    bool flag = false;

    if (C <= 0)
        return scalingFactor;
    else if (C >= 1)
        return 0;
    else if (C > 0.5)
    {
        C = 1 - C;
        flag = true;
    }

    mm = Foam::min(m12, m3);

    // Small modification to prevent round off errors in some extreme cases
    // where m1 is close to 0.

    if ((m2 <= C) & (C < mm) & ((m1/m2) < 1e-12))
    {
        m1 = 0;
        m12 = m2;
        mm = Foam::min(m12, m3);
    }


    scalar v1 = Foam::sqr(m1) / Foam::max(6 * m2 * m3, 1e-50);

    if (C < m1)
    {
        V = Foam::pow3(C) / (6 * m1 * m2 * m3);
    }
    else if (C < m2)
    {
        V = (C*(C - m1)) / (2*m2*m3) + v1;
    }
    else if (C < mm)
    {
        V = (Foam::sqr(C) * (3*m12 - C) + Foam::sqr(m1) * (m1 - 3*C) + Foam::sqr(m2) * (m2 - 3*C))
            / (6*m1*m2*m3);
    }
    else if (m3 < m12)
    {
        V = (Foam::sqr(C) * (3 - 2*C) + Foam::sqr(m1) * (m1 - 3*C)
            + Foam::sqr(m2) * (m2 - 3*C) + Foam::sqr(m3) * (m3 - 3*C))
            / (6*m1*m2*m3);
    }
    else
    {
        V = (2*C - m12) / (2*m3);
    }


    if (flag)
    {
        return (V) * scalingFactor;
    }
    else
    {
        return (1-V) * scalingFactor;
    }

}

}

}

}
