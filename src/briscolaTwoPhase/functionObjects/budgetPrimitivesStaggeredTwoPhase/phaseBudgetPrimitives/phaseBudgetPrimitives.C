#include "phaseBudgetPrimitives.H"
#include "budgetPrimitivesStaggeredTwoPhase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(phaseBudgetPrimitives, 0);

phaseBudgetPrimitives::phaseBudgetPrimitives
(
    const budgetPrimitivesStaggeredTwoPhase& bp,
    const colocatedScalarField& f,
    const label phase
)
:
    bp_(bp),
    fvMsh_(bp_.fvMsh_),
    rho_(bp_.tpm_.dict().lookup<scalar>("rho"+Foam::name(phase))),
    mu_(bp_.tpm_.dict().lookup<scalar>("mu"+Foam::name(phase))),
    f_(f),
    Uc_(bp_.Uc_),
    p_(bp_.p_),
    si_(bp_.si_),
    Sij_(bp_.Sij_),
    djui_(bp_.djui_),
    djf_
    (
        "djf"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    gamma_
    (
        "gamma"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    djG_
    (
        "djG"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    fui_
    (
        "fui"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    fuiui_
    (
        "fuiui"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    djfuiui_
    (
        "djfuiui"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    fuiuj_
    (
        "fuiuj"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    uidjf_
    (
        "uidjf"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    djfui_
    (
        "djfui"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    fdjui_
    (
        "fdjui"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    djfuiuj_
    (
        "djfuiuj"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    djfuiuiuj_
    (
        "djfuiuiuj"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    djfuiuiujalt_
    (
        "djfuiuiuj_alt"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    fQ_
    (
        "fQ"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    djfQ_
    (
        "djfQ"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    fQui_
    (
        "fQui"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    djfQui_
    (
        "djfQui"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    uidjfnuSij_
    (
        "uidjfnuSij"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    djfnuSxj_
    (
        "djfnuSxj"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    djfnuSyj_
    (
        "djfnuSyj"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    djfnuSzj_
    (
        "djfnuSzj"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    fnuSij_
    (
        "fnuSij"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    fnuSxj_
    (
        "fnuSxj"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    fnuSyj_
    (
        "fnuSyj"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    fnuSzj_
    (
        "fnuSzj"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    fnuSijdjui_
    (
        "fnuSijdjui"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    fTij_
    (
        "fTij"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    rhoTijdjG_
    (
        "rhoTijdjG"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    rhoTijuidjG_
    (
        "rhoTijuidjG"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    fsi_
    (
        "fsi"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    fuisi_
    (
        "fuisi"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    ff_
    (
        "ff"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    Tij_
    (
        "Tij"+Foam::name(phase),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    )
{
    djf_           = Zero;
    gamma_         = Zero;
    djG_           = Zero;
    fui_           = Zero;
    fuiuj_         = Zero;
    fuiui_         = Zero;
    djfuiui_       = Zero;
    uidjf_         = Zero;
    djfui_         = Zero;
    fdjui_         = Zero;
    djfuiuj_       = Zero;
    djfuiuiuj_     = Zero;
    djfuiuiujalt_  = Zero;
    fQ_            = Zero;
    djfQ_          = Zero;
    fQui_          = Zero;
    djfQui_        = Zero;
    uidjfnuSij_    = Zero;
    djfnuSxj_      = Zero;
    djfnuSyj_      = Zero;
    djfnuSzj_      = Zero;
    fnuSij_        = Zero;
    fnuSxj_        = Zero;
    fnuSyj_        = Zero;
    fnuSzj_        = Zero;
    fTij_          = Zero;
    rhoTijdjG_     = Zero;
    rhoTijuidjG_   = Zero;
    fsi_           = Zero;
    fuisi_         = Zero;
    ff_            = Zero;
    Tij_           = Zero;
}

phaseBudgetPrimitives::~phaseBudgetPrimitives()
{}

void phaseBudgetPrimitives::correct()
{
    const scalar nu = mu_ / rho_;

    gamma_ = f_ / rho_;
    gamma_.correctBoundaryConditions();

    djG_ = ex::grad(gamma_);

    const colocatedVectorField& CC =
        fvMsh_.metrics<colocated>().cellCenters();

    djf_ = ex::grad(f_);

    fui_ = f_*Uc_;

    djfui_ = ex::grad(fui_);

    fuiui_ = f_*(Uc_ & Uc_);

    djfuiui_ = ex::grad(fuiui_);

    fuiuj_ = f_*Uc_*Uc_;

    uidjf_ = Uc_*djf_;

    fdjui_ = f_ * ex::grad(Uc_);

    // Product rule fir djfuiuj:

    // forAllCells(djfuiuj_,i,j,k)
    // {
    //     const vector& Uc = Uc_(i,j,k);
    //     const tensor& djui = djui_(i,j,k);
    //     const vector& djf = djf_(i,j,k);

    //     vector& djfuiuj = djfuiuj_(i,j,k);

    //     for (int m = 0; m < 3; m++)
    //         djfuiuj[m] = Uc[m] * Uc.x() * djf.x()
    //                    + Uc[m] * Uc.y() * djf.y()
    //                    + Uc[m] * Uc.z() * djf.z()
    //                    + f_(i,j,k) *
    //                    (
    //                       Uc.x() * djui(m,0)
    //                     + Uc.y() * djui(m,1)
    //                     + Uc.z() * djui(m,2)
    //                    )
    //                    + f_(i,j,k) * Uc[m] *
    //                    (
    //                       djui.xx()
    //                     + djui.yy()
    //                     + djui.zz()
    //                    );
    // }

    // Divergence operator for djfuiuj:

    djfuiuj_ = ex::div(ex::faceFlux(fuiuj_));

    fQ_ = f_ * ((p_  / rho_) + (bp_.tpm_.g() & CC));

    djfQ_ = ex::grad(fQ_);

    fQui_ = fQ_ * Uc_;

    djfQui_ = ex::grad(fQ_*Uc_);

    forAllCells(Tij_,i,j,k)
    {
        scalar gh = bp_.tpm_.g() & CC(i,j,k);

        Tij_(i,j,k) = 2.0 * nu * Sij_(i,j,k)
            - tensor::I * (p_(i,j,k) / rho_ + gh);
    }

    djfuiuiuj_ =
          ex::div(ex::faceFlux(fuiuj_.component(tensor::XX)*Uc_))
        + ex::div(ex::faceFlux(fuiuj_.component(tensor::YY)*Uc_))
        + ex::div(ex::faceFlux(fuiuj_.component(tensor::ZZ)*Uc_));

    forAllCells(djfuiuiujalt_,i,j,k)
    {
        const vector& Uc = Uc_(i,j,k);
        const vector& djf = djf_(i,j,k);
        const tensor& djui = djui_(i,j,k);

        djfuiuiujalt_(i,j,k) = f_(i,j,k)
                * (Foam::sqr(Uc.x()) + Foam::sqr(Uc.y()) + Foam::sqr(Uc.z()))
                * (djui.xx() + djui.yy() + djui.zz())
                + 2.0 * f_(i,j,k) *
                (
                      Uc.x() * Uc.x() * djui.xx()
                    + Uc.x() * Uc.y() * djui.xy()
                    + Uc.x() * Uc.z() * djui.xz()
                    + Uc.y() * Uc.x() * djui.yx()
                    + Uc.y() * Uc.y() * djui.yy()
                    + Uc.y() * Uc.z() * djui.yz()
                    + Uc.z() * Uc.x() * djui.zx()
                    + Uc.z() * Uc.y() * djui.zy()
                    + Uc.z() * Uc.z() * djui.zz()
                )
                + Foam::pow3(Uc.x()) * djf.x()
                + Foam::pow3(Uc.y()) * djf.y()
                + Foam::pow3(Uc.z()) * djf.z();
    }

    fnuSij_ = f_ * nu * Sij_;

    fnuSxj_ = (T(fnuSij_) & vector(1,0,0));
    fnuSyj_ = (T(fnuSij_) & vector(0,1,0));
    fnuSzj_ = (T(fnuSij_) & vector(0,0,1));

    djfnuSxj_ = ex::div(ex::faceFlux(fnuSxj_));
    djfnuSyj_ = ex::div(ex::faceFlux(fnuSyj_));
    djfnuSzj_ = ex::div(ex::faceFlux(fnuSzj_));

    fnuSijdjui_ = (fnuSij_ && djui_);

    uidjfnuSij_ =
        Uc_.component(vector::X)*djfnuSxj_
      + Uc_.component(vector::Y)*djfnuSyj_
      + Uc_.component(vector::Z)*djfnuSzj_;

    fTij_ = f_ * Tij_;

    fsi_ = f_ * si_ / rho_;

    ff_ = f_ * f_;

    fuisi_ = (Uc_ & si_);

    rhoTijdjG_ = rho_ * (Tij_ & djG_);

    rhoTijuidjG_ = rho_ * (Tij_ && (Uc_*djG_));
}

}

}

}

}
