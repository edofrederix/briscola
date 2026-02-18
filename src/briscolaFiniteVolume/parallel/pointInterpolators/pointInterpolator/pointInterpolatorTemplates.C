#include "pointInterpolator.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Operator that sets x to y if x is zero

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

template<class MeshType>
template<class Type>
List<Type> pointInterpolator<MeshType>::combine(List<Type>& values) const
{
    Pstream::listCombineGather(values, nonZeroOp<Type>());
    Pstream::listCombineScatter(values);

    if (global_)
    {
        return List<Type>(values);
    }
    else
    {
        return List<Type>(SubList<Type>(values, size_, start_));
    }
}

}

}

}
