namespace Foam
{

namespace briscola
{

namespace fv
{

inline truncatedPiped::truncatedPiped
(
    const vector& lba,
    const vector& rba,
    const vector& lta,
    const vector& lbf,
    const vector& n,
    const scalar& C
)
:
    lba_(lba),
    rba_(rba),
    lta_(lta),
    lbf_(lbf),
    n_(n),
    C_(C)
{}

inline truncatedPiped::truncatedPiped
(
    const vertexVector& v,
    const vector& n,
    const scalar& C
)
:
    lba_(v.lba()),
    rba_(v.rba()),
    lta_(v.lta()),
    lbf_(v.lbf()),
    n_(n),
    C_(C)
{}

inline truncatedPiped::truncatedPiped(const truncatedPiped& p)
:
    lba_(p.lba_),
    rba_(p.rba_),
    lta_(p.lta_),
    lbf_(p.lbf_),
    n_(p.n_),
    C_(p.C_)
{}

inline scalar truncatedPiped::volume() const
{
    // This function computes and returns the volume of the parallelepiped
    // truncated by the plane dot(x,n) + C = 0.
    //
    // First the linear transform
    //
    // x_old = lba + T * x_new
    //
    // projects the parallelepiped into the unit cube, then the algorithm of
    // Scardovelli & Zaleski (2000) is used.

    const tensor T
    (
        rba_ - lba_,
        lta_ - lba_,
        lbf_ - lba_
    );

    vector n = T & n_;
    scalar C = C_ + (n_ & lba_);

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

    if (C <= 0.0)
    {
        return scalingFactor;
    }
    else if (C >= 1.0)
    {
        return 0;
    }

    // Solve inverse problem for C > 0.5

    bool inverse = false;

    if (C > 0.5)
    {
        C = 1.0 - C;
        inverse = true;
    }

    mm = Foam::min(m12, m3);

    // Small modification to prevent round off errors in some extreme cases
    // where m1 is close to 0.

    if ((m2 <= C) && (C < mm) && ((m1/m2) < 1e-12))
    {
        m1 = 0;
        m12 = m2;
        mm = Foam::min(m12, m3);
    }

    scalar v1 = Foam::sqr(m1) / Foam::max(6.0 * m2 * m3, 1e-50);

    if (C < m1)
    {
        V = Foam::pow3(C) / (6.0 * m1 * m2 * m3);
    }
    else if (C < m2)
    {
        V = (C*(C - m1)) / (2.0*m2*m3) + v1;
    }
    else if (C < mm)
    {
        V =
            (
                Foam::sqr(C) * (3.0*m12 - C)
              + Foam::sqr(m1) * (m1 - 3.0*C)
              + Foam::sqr(m2) * (m2 - 3.0*C)
            )
          / (6*m1*m2*m3);
    }
    else if (m3 < m12)
    {
        V =
            (
                Foam::sqr(C)  * (3.0 - 2.0*C)
              + Foam::sqr(m1) * (m1 - 3.0*C)
              + Foam::sqr(m2) * (m2 - 3.0*C)
              + Foam::sqr(m3) * (m3 - 3.0*C)
            )
          / (6*m1*m2*m3);
    }
    else
    {
        V = (2.0*C - m12) / (2.0*m3);
    }

    if (inverse)
    {
        return V * scalingFactor;
    }
    else
    {
        return (1.0-V) * scalingFactor;
    }
}

}

}

}
