#include "LVE.H"
#include "vof.H"
#include "truncatedHex.H"
#include "SortableList.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

const scalar LVE::tol_ = 1e-8;
const label LVE::maxIter_ = 100;

scalar LVE::solveBracketedLVE
(
    const vertexVector& v,
    const scalar& V,
    const vector& n,
    const scalar& f,
    scalar CMin,
    scalar CMax,
    scalar VMin,
    scalar VMax
) const
{
    // Here we should exploit the analytical solution method of Lopez &
    // Hernandez (JCP, 2008). However, for now we've implemented a simple
    // bisection method assuming linear behavior between CMin and CMax.

    label iter = 0;

    const scalar Vf = f*V;

    scalar Vi = -V;
    scalar C = (CMin+CMax)/2.0;
    scalar pos = 0.5;

    while (Foam::mag((Vf-Vi)/V) > tol_ && iter < maxIter_)
    {
        // Use bisection of we are close to the edges of the [VMin,VMax]
        // interval. Use a linear interpolation otherwise.

        if (pos < 1e-3 || pos > 0.999)
        {
            C = (CMin+CMax)/2.0;
        }
        else
        {
            C = CMin + (Vf-VMin)/(VMax-VMin) * (CMax-CMin);
        }

        Vi = truncatedHex(v,n,C).volume();
        pos = (Vi - VMin)/(VMax-VMin+1e-50);

        #ifdef FULLDEBUG

        if (Vi - VMin < -1e-12 || VMax - Vi < -1e-12)
        {
            FatalErrorInFunction
                << "Bracketed LVE problem solver failed because the "
                << "truncated hex volume is non-monotonic." << endl
                << abort(FatalError);
        }

        #endif

        // Adjust search boundaries

        if (Vi < Vf)
        {
            CMin = C;
            VMin = Vi;
        }
        else
        {
            CMax = C;
            VMax = Vi;
        }

        iter++;
    }

    #ifdef FULLDEBUG

    if (iter == maxIter_)
    {
        WarningInFunction
            << "Bracketed LVE problem solver did not converge" << endl;
    }

    #endif

    return C;
}

LVE::LVE(const vof& vf)
:
    vf_(vf),
    fvMsh_(vf.fvMsh()),
    rectilinear_(fvMsh_.msh()[0].rectilinear() == unitXYZ)
{}

LVE::~LVE()
{}

scalar LVE::operator()(const labelVector& ijk, const vector& n) const
{
    const scalar f(vf_.alpha()[0][0](ijk));

    const vertexVector vertices
    (
        fvMsh_.template metrics<colocated>().vertexCenters()[0][0](ijk)
    );

    const vector m(n/Foam::mag(n));

    return rectilinear_ ? pLVE(vertices,m,f) : aLVE(vertices,m,f);
}

scalar LVE::aLVE
(
    const vertexVector& v,
    const vector& n,
    const scalar& f
) const
{
    scalarList C(8);
    scalarList V(8);

    for (int i = 0; i < 8; i++)
    {
        C[i] = -(n & v[i]);
        V[i] = truncatedHex(v,n,C[i]).volume();
    }

    labelList order;
    uniqueOrder(V, order);

    const scalar Vc = V[order[order.size()-1]];

    // Early detection at the vertices

    forAll(order, i)
    {
        if (Foam::mag(f-V[order[i]]/Vc) < tol_)
        {
            return C[order[i]];
        }
    }

    const scalar Vf = Vc*f;

    // Determine bracket

    scalar CMin = C[order[0]];
    scalar CMax = C[order[order.size()-1]];
    scalar VMin = V[order[0]];
    scalar VMax = V[order[order.size()-1]];

    for (int i = 1; i < order.size(); i++)
    {
        if (V[order[i-1]] < Vf && Vf < V[order[i]])
        {
            CMin = C[order[i-1]];
            CMax = C[order[i]];

            VMin = V[order[i-1]];
            VMax = V[order[i]];

            break;
        }
    }

    // Solve

    return solveBracketedLVE(v,Vc,n,f,CMin,CMax,VMin,VMax);
}

scalar LVE::pLVE
(
    const vertexVector& vertices,
    const vector& no,
    const scalar& fo
) const
{
    // Algorithm of Scardovelli & Zaleski (2000) as implemented in Basilisk and
    // VOFTools

    #ifdef FULLDEBUG

    vector d
    (
        Foam::mag(vertices.rba() - vertices.lba()),
        Foam::mag(vertices.lta() - vertices.lba()),
        Foam::mag(vertices.lbf() - vertices.lba())
    );

    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
    {
        const label v0 = vertexNumsInEdge[i*4+j][0];
        const label v1 = vertexNumsInEdge[i*4+j][1];
        const scalar dj = Foam::mag(vertices[v0] - vertices[v1]);

        if (Foam::mag(dj - d[i]) > 1e-12)
            FatalErrorInFunction
                << "Cell is not a parallelepiped" << endl
                << abort(FatalError);
    }

    #endif

    const tensor T
    (
        vertices.rba() - vertices.lba(),
        vertices.lta() - vertices.lba(),
        vertices.lbf() - vertices.lba()
    );

    vector n = T & no;
    const scalar S = Foam::mag(n);
    n = n/S;

    scalarList Cs(8);

    for (int i = 0; i < 8; i++)
    {
        Cs[i] = -(n & unitCube[i]);
    }

    const scalar CMin = Foam::min(Cs);
    const scalar CMax = Foam::max(Cs);

    // Solve inverse if f > 0.5

    const scalar f = Foam::max(Foam::min(fo, 1.0-fo),0);

    scalar s = cmptSum(cmptMag(n));
    vector m = cmptMag(n)/s;

    scalarList ml(3);

    for (int i = 0; i < 3; i++)
        ml[i] = m[i];

    sort(ml);

    const scalar& m1 = ml[0];
    const scalar& m2 = ml[1];
    const scalar& m3 = ml[2];

    const scalar m12 = m1 + m2;

    const scalar pr = Foam::max(6.0*m1*m2*m3, 1e-50);
    const scalar v1 = Foam::pow3(m1)/pr;
    const scalar v2 = v1 + (m2 - m1)/(2.0*m3);

    scalar v3, mm;

    if (m3 < m12)
    {
        mm = m3;
        v3 =
            (
              + (3*m12 - m3)*Foam::sqr(m3)
              + (m1 - 3*m3 )*Foam::sqr(m1)
              + (m2 - 3*m3 )*Foam::sqr(m2)
            )
          / pr;
    }
    else
    {
        mm = m12;
        v3 = m12/(2.0*m3);
    }

    scalar alpha = 0;

    if (f < v1)
    {
        alpha = Foam::pow(pr*f, 1.0/3.0);
    }
    else if (f < v2)
    {
        alpha = 0.5 * (m1 + Foam::sqrt(Foam::sqr(m1) + 8.0*m2*m3*(f - v1)));
    }
    else if (f < v3)
    {
        const scalar p12 = Foam::sqrt(2.0*m1*m2);
        const scalar q = 3.0*(m12 - 2.0*m3*f)/(4.0*p12);
        const scalar theta = Foam::acos(Foam::max(Foam::min(q,1),-1))/3.0;
        const scalar cs = Foam::cos(theta);

        alpha = p12*(Foam::sqrt(3.0*(1.0 - Foam::sqr(cs))) - cs) + m12;
    }
    else if (m12 <= m3)
    {
        alpha = m3*f + 0.5*mm;
    }
    else
    {
        const scalar p = m1*(m2 + m3) + m2*m3 - 0.25;
        const scalar p12 = Foam::sqrt(p);
        const scalar q = 3.0*m1*m2*m3*(0.5 - f)/(2.0*p*p12);
        const scalar theta = Foam::acos(Foam::max(Foam::min(q,1),-1))/3.0;
        const scalar cs = Foam::cos(theta);

        alpha = p12*(Foam::sqrt(3.0*(1.0 - Foam::sqr(cs))) - cs) + 0.5;
    }

    scalar C;

    if (fo <= 0.5)
    {
        C = CMin + alpha * Foam::mag(CMax - CMin);
    }
    else
    {
        C = CMax - alpha * Foam::mag(CMax - CMin);
    }

    return C*S - (no & vertices.lba());
}

}

}

}
