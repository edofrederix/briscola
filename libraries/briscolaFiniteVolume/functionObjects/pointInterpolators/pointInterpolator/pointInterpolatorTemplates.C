#include "pointInterpolator.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class nonZeroOp
{
    public:

        void operator()(Type& x, Type& y) const
        {
            if (x == pTraits<Type>::zero)
            {
                x = y;
            }
        }
};

template<class Type>
void pointInterpolator::gatherScatter(List<Type>& values)
{
    Pstream::listCombineGather(values, nonZeroOp<Type>());
    Pstream::listCombineScatter(values);
}

}

}

}
