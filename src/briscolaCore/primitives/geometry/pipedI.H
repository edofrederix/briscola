namespace Foam
{

namespace briscola
{

inline scalar piped::truncationVolume(const vector& m, const scalar C) const
{
    // This function computes the volume of the piped truncated by the plane
    // dot(x,n) + C = 0 using the algorithm of Scardovelli & Zaleski (2000).
    // This is faster than decomposing the piped in six tets and truncating
    // those.

    vector n = T() & m;
    scalar D = C + (m & p());

    for (label i = 0; i < 3; i++)
    {
        if (n[i] < 0)
        {
            D += n[i];
            n[i] = -n[i];
        }
    }

    const scalar f = volume();
    const scalar s = cmptSum(cmptMag(n));

    n = cmptMag(n)/s;
    D = -D/s;

    scalarList ml(3);

    for (int i = 0; i < 3; i++)
        ml[i] = n[i];

    sort(ml);

    scalar& m1 = ml[0], m2 = ml[1], m3 = ml[2];

    scalar m12 = m1 + m2;
    scalar V, mm;

    if (D <= 0.0)
        return f;

    if (D >= 1.0)
        return 0.0;

    // Solve inverse problem for D > 0.5

    bool inverse = false;

    if (D > 0.5)
    {
        D = 1.0 - D;
        inverse = true;
    }

    mm = Foam::min(m12, m3);

    // Small modification to prevent round off errors in some extreme cases
    // where m1 is close to 0.

    if (m2 <= D && D < mm && m1/m2 < 1e-12)
    {
        m1 = 0;
        m12 = m2;
        mm = Foam::min(m12, m3);
    }

    scalar v1 = Foam::sqr(m1)/Foam::max(6.0*m2*m3, 1e-50);

    if (D < m1)
    {
        V = Foam::pow3(D)/(6.0*m1*m2*m3);
    }
    else if (D < m2)
    {
        V = D*(D - m1)/(2.0*m2*m3) + v1;
    }
    else if (D < mm)
    {
        V =
            (
                Foam::sqr(D)*(3.0*m12 - D)
              + Foam::sqr(m1)*(m1 - 3.0*D)
              + Foam::sqr(m2)*(m2 - 3.0*D)
            )
        / (6.0*m1*m2*m3);
    }
    else if (m3 < m12)
    {
        V =
            (
                Foam::sqr(D)*(3.0 - 2.0*D)
              + Foam::sqr(m1)*(m1 - 3.0*D)
              + Foam::sqr(m2)*(m2 - 3.0*D)
              + Foam::sqr(m3)*(m3 - 3.0*D)
            )
        / (6.0*m1*m2*m3);
    }
    else
    {
        V = 0.5*(2.0*D - m12)/m3;
    }

    if (inverse)
    {
        return V*f;
    }
    else
    {
        return (1.0 - V)*f;
    }
}

}

}
