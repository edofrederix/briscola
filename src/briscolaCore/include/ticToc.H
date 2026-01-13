// Some convenient macro's that help with profiling

#include <chrono>

#define initTicTocHeader scalarList ticTocs;
#define initTicTocConstructor(count)  ticTocs(count,Zero)

#define initTicToc(count) scalarList ticTocs(count,Zero);

#define tic(num)                                                               \
    auto ticToc##num = ::std::chrono::high_resolution_clock::now();

#define toc(num)                                                               \
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
    Info<< "Tic-tocs [ms] = " << endl;                                         \
    forAll(ticTocs, i)                                                         \
        Info<< Foam::name(i) << "=" << ticTocs[i]/1000.0 << " ("               \
            << ticTocs[i]/sum(ticTocs)*100 << "%)" << endl;                    \
    Info<< "sum = " << sum(ticTocs)/1000.0 << endl;                            \
}
