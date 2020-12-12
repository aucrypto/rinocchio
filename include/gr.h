#ifndef GR_H
#define GR_H

#include <NTL/ZZ_pE.h>

using namespace NTL;

// Random element in R*
ZZ_pE randomInvertible();

ZZ_pE indexedElementInExceptionalSet(long index);

Vec<ZZ_pE> getExceptionalSubset(long size);

// Random element in A
ZZ_pE randomInExceptionalSet();

// Random element in A*
ZZ_pE randomNonZeroInExceptionalSet();

ZZ_pE getInverse(ZZ_pE element);

ZZ_pX primitiveIrredPoly(long degree);

#endif