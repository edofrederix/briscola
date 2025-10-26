#ifndef linearSystemFunctions_H
#define linearSystemFunctions_H

#include "stencil.H"
#include "diagStencil.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace linearSystemFun
{

// Low level row product

    // Forward declaration

    template<class SType, class Type>
    class rowProduct;

    // Stencil

    template<class Type>
    class rowProduct<stencil,Type>
    {
    public:
        inline Type operator()
        (
            const scalar* const __restrict__ st_arr,
            const Type* const __restrict__ s_arr,
            const label i,
            const label S_i,
            const label S_j,
            const label S_k
        )
        {
            return
                st_arr[i*7  ] * s_arr[i      ]
              + st_arr[i*7+1] * s_arr[i - S_i]
              + st_arr[i*7+2] * s_arr[i + S_i]
              + st_arr[i*7+3] * s_arr[i - S_j]
              + st_arr[i*7+4] * s_arr[i + S_j]
              + st_arr[i*7+5] * s_arr[i - S_k]
              + st_arr[i*7+6] * s_arr[i + S_k];
        }
    };

    // Diagonal stencil

    template<class Type>
    class rowProduct<diagStencil,Type>
    {
    public:
        inline Type operator()
        (
            const scalar* const __restrict__ st_arr,
            const Type* const __restrict__ s_arr,
            const label i,
            const label S_i,
            const label S_j,
            const label S_k
        )
        {
            return st_arr[i]*s_arr[i];
        }
    };

// Low level off-diagonal row product

    // Forward declaration

    template<class SType, class Type>
    class offDiagRowProduct;

    // Stencil

    template<class Type>
    class offDiagRowProduct<stencil,Type>
    {
    public:
        inline Type operator()
        (
            const scalar* const __restrict__ st_arr,
            const Type* const __restrict__ s_arr,
            const label i,
            const label S_i,
            const label S_j,
            const label S_k
        )
        {
            return
                st_arr[i*7+1] * s_arr[i - S_i]
              + st_arr[i*7+2] * s_arr[i + S_i]
              + st_arr[i*7+3] * s_arr[i - S_j]
              + st_arr[i*7+4] * s_arr[i + S_j]
              + st_arr[i*7+5] * s_arr[i - S_k]
              + st_arr[i*7+6] * s_arr[i + S_k];
        }
    };

    // Diagonal stencil

    template<class Type>
    class offDiagRowProduct<diagStencil,Type>
    {
    public:
        inline Type operator()
        (
            const scalar* const __restrict__ st_arr,
            const Type* const __restrict__ s_arr,
            const label i,
            const label S_i,
            const label S_j,
            const label S_k
        )
        {
            return Zero;
        }
    };

}

}

}

}

#endif
