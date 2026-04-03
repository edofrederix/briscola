// Some convenient macro's that help with profiling

#include <chrono>

#define initTicTocHeader scalarList ticTocs;
#define initTicTocConstructor(count)  ticTocs(count,Zero)

#define initTicToc(count) scalarList ticTocs(count,Zero);

#define Tic(num)                                                               \
    auto ticToc##num = ::std::chrono::high_resolution_clock::now();

#define Toc(num)                                                               \
    ticTocs[num] +=                                                            \
        scalar                                                                 \
        (                                                                      \
            ::std::chrono::duration_cast<::std::chrono::microseconds>          \
            (                                                                  \
                ::std::chrono::high_resolution_clock::now()                    \
              - ticToc##num                                                    \
             ).count()                                                         \
        );

#define printTicToc                                                            \
{                                                                              \
    Info<< "Tic-toc parallel averages = " << endl;                             \
    forAll(ticTocs, i)                                                         \
        Info<< Foam::name(i) << "="                                            \
            << returnReduce(ticTocs[i],sumOp<scalar>())                        \
             / Pstream::nProcs()/1000.0 << " ms ("                             \
            << returnReduce(ticTocs[i],sumOp<scalar>())                        \
             / sum(returnReduce(ticTocs, ListOp<sumOp<scalar>>()))*100         \
            << "%)" << endl;                                                   \
    Info<< "Mean sum = "                                                       \
        << sum(returnReduce(ticTocs, ListOp<sumOp<scalar>>()))                 \
         / Pstream::nProcs()/1000.0 << " ms" << endl;                          \
}
