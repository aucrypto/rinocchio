#ifndef GR_H
#define GR_H

#include <NTL/ZZ_pE.h>

using namespace NTL;

// Random element in R*
ZZ_pE randomInvertible();

ZZ_pE indexedElementInExceptionalSet(long index);

// Random element in A
ZZ_pE randomInExceptionalSet();

// Random element in A*
ZZ_pE randomNonZeroInExceptionalSet();

#endif