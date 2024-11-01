#ifndef RungeKuttaSchemes_H
#define RungeKuttaSchemes_H

#include "RungeKuttaScheme.H"

#define makeRungeKuttaScheme(SCHEME)                                           \
    defineTypeNameAndDebug(SCHEME, 0)                                          \
    addToRunTimeSelectionTable(RungeKuttaScheme, SCHEME, dictionary)

#endif
