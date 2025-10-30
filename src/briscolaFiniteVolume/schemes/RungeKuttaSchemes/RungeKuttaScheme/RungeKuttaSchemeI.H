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
    return a().m();
}

inline scalar RungeKuttaScheme::A() const
{
    return a()(stage_-1,stage_-1);
}

inline scalar RungeKuttaScheme::B() const
{
    return b()(stage_-1,stage_-1);
}

inline scalar RungeKuttaScheme::C() const
{
    scalar sum = 0.0;
    for (int j = 0; j < a().n(); j++)
        sum += a()(stage_-1,j);

    return sum;
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
    for (int i = 0; i < a().m(); i++)
        if (a()(i,i))
            return false;

    return true;
}

inline Switch RungeKuttaScheme::exB() const
{
    for (int i = 0; i < b().m(); i++)
        if (b()(i,i))
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

inline Switch RungeKuttaScheme::storeStageA() const
{
    for (label i = stage_+1; i <= numberOfStages(); i++)
        if (a()(i-1,stage_-1))
            return true;

    return false;
}

inline Switch RungeKuttaScheme::storeStageB() const
{
    for (label i = stage_+1; i <= numberOfStages(); i++)
        if (b()(i-1,stage_-1))
            return true;

    return false;
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
