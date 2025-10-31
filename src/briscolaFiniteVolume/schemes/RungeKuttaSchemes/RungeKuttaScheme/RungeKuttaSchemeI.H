namespace Foam
{

namespace briscola
{

namespace fv
{

template<>
inline FastPtrList<FastPtrList<colocatedScalarField>>&
RungeKuttaScheme::stageSources<colocatedScalarField>()
{
    return coloScalarStageSources_;
}

template<>
inline FastPtrList<FastPtrList<colocatedVectorField>>&
RungeKuttaScheme::stageSources<colocatedVectorField>()
{
    return coloVectorStageSources_;
}

template<>
inline FastPtrList<FastPtrList<staggeredScalarField>>&
RungeKuttaScheme::stageSources<staggeredScalarField>()
{
    return stagScalarStageSources_;
}

template<>
inline FastPtrList<FastPtrList<staggeredVectorField>>&
RungeKuttaScheme::stageSources<staggeredVectorField>()
{
    return stagVectorStageSources_;
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
    return stage_ == nStages_;
}

inline Switch RungeKuttaScheme::storeStageA() const
{
    for (int step = 0; step < nSteps_; step++)
        for (label i = stage_; i < nStages_; i++)
            if (a()(i,stage_-1+step*nStages_))
                return true;

    return false;
}

inline Switch RungeKuttaScheme::storeStageB() const
{
    for (int step = 0; step < nSteps_; step++)
        for (label i = stage_; i < nStages_; i++)
            if (b()(i,stage_-1+step*nStages_))
                return true;

    return false;
}

inline Switch RungeKuttaScheme::solve() const
{
    return !(firstStage() && exStage());
}

template<class FieldType>
void RungeKuttaScheme::cycle()
{
    FastPtrList<FastPtrList<FieldType>>& sources =
        stageSources<FieldType>();

    // Move pointers forward by nStages and recycle the oldest step

    labelList oldToNew(nStages_*nSteps_);
    forAll(oldToNew, i)
        oldToNew[i] = (i + nStages_) % (nStages_*nSteps_);

    forAll(sources, i)
    {
        PtrList<FieldType>& list = sources[i];

        // Reorder

        list.reorder(oldToNew);

        // Rename and reset write option

        int j = 0;

        for (int step = 0; step < nSteps_; step++)
        {
            for (int stage = 0; stage < nStages_; stage++)
            {
                const word baseName = IOobject::member(list[j].name());

                if (step == 0)
                {
                    list[j].rename(IOobject::groupName(baseName, "0"));
                }
                else
                {
                    list[j].rename(list[j-nStages_].name() + "_0");
                }

                // Write everything but the last step

                if (step < nSteps_-1)
                {
                    list[j].writeOpt() = IOobject::AUTO_WRITE;
                }
                else
                {
                    list[j].writeOpt() = IOobject::NO_WRITE;
                }

                j++;
            }
        }
    }
}

inline Switch RungeKuttaScheme::loop()
{
    stage_++;

    if (stage_ > nStages_)
    {
        stage_ = 0;

        return false;
    }
    else
    {
        if (stage_ == 1 && nSteps_ > 1)
        {
            // For multi-step methods, cycle the sources

            cycle<colocatedScalarField>();
            cycle<colocatedVectorField>();
            cycle<staggeredScalarField>();
            cycle<staggeredVectorField>();
        }

        return true;
    }
}

}

}

}
