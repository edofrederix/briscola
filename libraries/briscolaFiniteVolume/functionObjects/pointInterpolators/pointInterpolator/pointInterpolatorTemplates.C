#include "pointInterpolator.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class averageNonZeroOp
{
    public:

        void operator()(Type& x, Type& y) const
        {
            if (x == pTraits<Type>::zero)
            {
                x = y;
            }
            else if (y != pTraits<Type>::zero)
            {
                x = (x+y)/2.0;
            }
        }
};

template<class Type>
void pointInterpolator::gatherScatter(List<Type>& values)
{
    Pstream::listCombineGather(values, averageNonZeroOp<Type>());
    Pstream::listCombineScatter(values);
}

}

}

}
