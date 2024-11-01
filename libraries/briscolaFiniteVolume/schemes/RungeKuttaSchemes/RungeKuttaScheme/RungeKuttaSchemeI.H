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

inline Switch RungeKuttaScheme::explicitStageA() const
{
    return A() == 0.0;
}

inline Switch RungeKuttaScheme::explicitStageB() const
{
    return B() == 0.0;
}

inline Switch RungeKuttaScheme::explicitStage() const
{
    return explicitStageA() && explicitStageB();
}

inline Switch RungeKuttaScheme::implicitStageA() const
{
    return A() != 0.0;
}

inline Switch RungeKuttaScheme::implicitStageB() const
{
    return B() != 0.0;
}

inline Switch RungeKuttaScheme::implicitStage() const
{
    return implicitStageA() || implicitStageB();
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
    return !(firstStage() && explicitStage());
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
