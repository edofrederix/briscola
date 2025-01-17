namespace Foam
{

namespace briscola
{

namespace fv
{

inline const fvMesh& RungeKuttaScheme::fvMsh() const
{
    return fvMsh_;
}

inline const dictionary& RungeKuttaScheme::dict() const
{
    return dict_;
}

inline label RungeKuttaScheme::stage() const
{
    return stage_;
}

inline label RungeKuttaScheme::numberOfStages() const
{
    return a().size();
}

inline scalar RungeKuttaScheme::A() const
{
    return a()[stage_-1][stage_-1];
}

inline scalar RungeKuttaScheme::B() const
{
    return b()[stage_-1][stage_-1];
}

inline scalar RungeKuttaScheme::C() const
{
    return sum(a()[stage_-1]);
}

inline Switch RungeKuttaScheme::exStageA() const
{
    return A() == 0.0;
}

inline Switch RungeKuttaScheme::exStageB() const
{
    return B() == 0.0;
}

inline Switch RungeKuttaScheme::exStage() const
{
    return exStageA() && exStageB();
}

inline Switch RungeKuttaScheme::exA() const
{
    forAll(a(), i)
        if (a()[i][i])
            return false;

    return true;
}

inline Switch RungeKuttaScheme::exB() const
{
    forAll(b(), i)
        if (b()[i][i])
            return false;

    return true;
}

inline Switch RungeKuttaScheme::ex() const
{
    return exA() && exB();
}

inline Switch RungeKuttaScheme::imStageA() const
{
    return !exStageA();
}

inline Switch RungeKuttaScheme::imStageB() const
{
    return !exStageB();
}

inline Switch RungeKuttaScheme::imStage() const
{
    return imStageA() || imStageB();
}

inline Switch RungeKuttaScheme::imA() const
{
    return !exA();
}

inline Switch RungeKuttaScheme::imB() const
{
    return !exB();
}

inline Switch RungeKuttaScheme::im() const
{
    return imA() && imB();
}

inline Switch RungeKuttaScheme::firstStage() const
{
    return stage_ == 1;
}

inline Switch RungeKuttaScheme::lastStage() const
{
    return stage_ == numberOfStages();
}

inline Switch RungeKuttaScheme::solve() const
{
    return !(firstStage() && exStage());
}

inline Switch RungeKuttaScheme::loop()
{
    stage_++;

    if (stage_ > numberOfStages())
    {
        stage_ = 0;

        return false;
    }
    else
    {
        return true;
    }
}

}

}

}
