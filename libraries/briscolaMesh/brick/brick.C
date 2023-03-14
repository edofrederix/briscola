#include "brick.H"
#include "geometry.H"
#include "block.H"

namespace Foam
{

namespace briscola
{

void brick::createFaces()
{
    faces_.clear();
    faces_.setSize(6);

    // Walk through faces in a way that matches the faceNumsInBlock array

    label l = 0;

    for (label dir = 0; dir < 3; dir++)
    {
        for (label i = 0; i < 2; i++)
        {
            labelVector N(N_);
            N[dir] = 1;

            faces_.set
            (
                l,
                new face
                (
                    *this,
                    l,
                    v_.slice(i,dir,true)(),
                    N
                )
            );

            l++;
        }
    }
}

brick::brick(const geometry& g, const label num, const dictionary& dict)
:
    meshObject<geometry>(g, num),
    dict_(dict),
    v_(dict.lookup("vertices")),
    N_(dict.lookup("N")),
    grading_(grading::New(*this)),
    faces_()
{
    for (label i = 0; i < v_.size()-1; i++)
    {
        for (label j = i+1; j < v_.size(); j++)
        {
            if (v_(i) == v_(j))
            {
                FatalErrorInFunction
                    << "Duplicate vertices found for " << *this
                    << exit(FatalError);
            }
        }
    }

    createFaces();

    if (leftHanded())
        FatalErrorInFunction
            << *this << " is left-handed. Bricks should be right-handed."
            << exit(FatalError);
}

brick::brick(const brick& b)
:
    meshObject<geometry>(b.parentGeometry(), b.num()),
    dict_(b.dict_),
    v_(b.v_),
    N_(b.N_),
    grading_(grading::New(*this)),
    faces_(b.faces_, *this)
{}

brick::brick(const brick& b, const geometry& geo)
:
    meshObject<geometry>(geo, b.num()),
    dict_(b.dict_),
    v_(b.v_),
    N_(b.N_),
    grading_(grading::New(*this)),
    faces_(b.faces_, *this)
{}

brick::~brick()
{}

void brick::transform(const labelTensor T)
{
    if (T == eye)
    {
        return;
    }
    else
    {
        v_.transform(T);

        N_ = cmptMag(T & N_);

        createFaces();

        grading_->transform(T);
    }
}

labelTensor brick::rotate(const label st, const label ax)
{
    const label s = (4+st%4)%4;

    labelTensor T = rotations[ax][s];

    this->transform(T);

    return T;
}

tmp<vectorBlock> brick::TFI
(
    const scalarField& xi,
    const scalarField& eta,
    const scalarField& zeta
) const
{
    const grading& g = grading_();

    // Gradings along 12 edges

    const scalarField g0(g(xi,0));
    const scalarField g1(g(xi,1));
    const scalarField g2(g(xi,2));
    const scalarField g3(g(xi,3));

    const scalarField g4(g(eta,4));
    const scalarField g5(g(eta,5));
    const scalarField g6(g(eta,6));
    const scalarField g7(g(eta,7));

    const scalarField g8(g(zeta,8));
    const scalarField g9(g(zeta,9));
    const scalarField g10(g(zeta,10));
    const scalarField g11(g(zeta,11));

    // Points along 12 edges given graded coordinates

    const vectorField e0g(e0()(g0));
    const vectorField e1g(e1()(g1));
    const vectorField e2g(e2()(g2));
    const vectorField e3g(e3()(g3));

    const vectorField e4g(e4()(g4));
    const vectorField e5g(e5()(g5));
    const vectorField e6g(e6()(g6));
    const vectorField e7g(e7()(g7));

    const vectorField e8g(e8()(g8));
    const vectorField e9g(e9()(g9));
    const vectorField e10g(e10()(g10));
    const vectorField e11g(e11()(g11));

    tmp<vectorBlock> tp(new vectorBlock(xi.size(), eta.size(), zeta.size()));

    vectorBlock& p = tp.ref();

    forAllBlock(p, i, j, k)
    {
        // Points on the graded unit cube

        const scalar xig =
            ((1.0-eta[j])*g0[i] + eta[j]*g1[i])*(1.0-zeta[k])
          + ((1.0-eta[j])*g2[i] + eta[j]*g3[i])*zeta[k];

        const scalar etag =
            ((1.0-zeta[k])*g4[j] + zeta[k]*g6[j])*(1.0-xi[i])
          + ((1.0-zeta[k])*g5[j] + zeta[k]*g7[j])*xi[i];

        const scalar zetag =
            ((1.0-xi[i])*g8[k] + xi[i]*g9[k])*(1.0-eta[j])
          + ((1.0-xi[i])*g10[k] + xi[i]*g11[k])*eta[j];

        // Transfinite interpolation

        p(i,j,k) =

            (1.0-etag) * (1.0-xig)   * e8g[k]
          + (1.0-etag) * xig         * e9g[k]
          + etag       * (1.0-xig)   * e10g[k]
          + etag       * xig         * e11g[k]

          + (1.0-etag) * (1.0-zetag) * e0g[i]
          + etag       * (1.0-zetag) * e1g[i]
          + (1.0-etag) * zetag       * e2g[i]
          + etag       * zetag       * e3g[i]

          + (1.0-xig) * (1.0-zetag)  * e4g[j]
          + xig       * (1.0-zetag)  * e5g[j]
          + (1.0-xig) * zetag        * e6g[j]
          + xig       * zetag        * e7g[j]

          - 2.0
          * (
                (1.0-xig) * (1.0-etag) * (1.0-zetag) * v0()
              + xig       * (1.0-etag) * (1.0-zetag) * v1()
              + (1.0-xig) * etag       * (1.0-zetag) * v2()
              + xig       * etag       * (1.0-zetag) * v3()
              + (1.0-xig) * (1.0-etag) * zetag       * v4()
              + xig       * (1.0-etag) * zetag       * v5()
              + (1.0-xig) * etag       * zetag       * v6()
              + xig       * etag       * zetag       * v7()
            );
    }

    return tp;
}

bool brick::operator==(const brick& b) const
{
    for (label i = 0; i < 8; i++)
    {
        for (label j = 0; j < 8; j++)
        {
            if (v_(i) == b.v()(j))
            {
                break;
            }
            else if (j == 7)
            {
                return false;
            }
        }
    }

    return true;
}

bool brick::operator!=(const brick& b) const
{
    return !(*this == b);
}

Ostream& operator<<(Ostream& os, const brick& b)
{
    os  << "brick " << b.num() << ": "
        << "(((" << b.v()(0) << "," << b.v()(1) << "),"
        << "(" << b.v()(2) << "," << b.v()(3) << "))"
        << "((" << b.v()(4) << "," << b.v()(5) << "),"
        << "(" << b.v()(6) << "," << b.v()(7) << ")))";

    return os;
}

}

}
