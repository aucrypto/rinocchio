#include <NTL/ZZ_pE.h>

using namespace NTL;

// Random element in R*
ZZ_pE randomInvertible() {
    while (true) {
        // random element in R:
        ZZ_pE res = random_ZZ_pE();
        //Check at least one coefficient is one
        for (int i = 0; i < ZZ_pE::degree(); i++) {
            //todo When d sufficiently large this check is not needed
            if (IsOdd(rep(rep(res)[i]))) return res;
            
        }
    }
}

ZZ_pE indexedElementInExceptionalSet(long index) {
    //e.g. 5=b101 becomes [1 0 1] and 11 = b1011 becomes [1 1 0 1]
    ZZ_pX res = ZZ_pX();
    for (int i = 0; i < ZZ_pE::degree(); i++) {
        long mask = 1 << i;
        if ((mask & index) != 0) {
            SetCoeff(res, i, 1);
        }
    }
    return to_ZZ_pE(res);
}

// Random element in A
ZZ_pE randomInExceptionalSet() {
    ZZ_pX a = ZZ_pX();
    for (int i = 0; i < ZZ_pE::degree(); i++) {
        long coeff;
        RandomBits(coeff, 1);
        if (coeff == 1) {
            SetCoeff(a, i, coeff);
        }
    }
    ZZ_pE fromPX = to_ZZ_pE(a);
    return fromPX;
}

// Random element in A*
ZZ_pE randomNonZeroInExceptionalSet() {
    while (true) {
        ZZ_pE res = randomInExceptionalSet();
        if(! IsZero(res)) return res;
    }
}