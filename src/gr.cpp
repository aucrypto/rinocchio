#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZX.h>

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
            SetCoeff(res, i);
        }
    }
    return to_ZZ_pE(res);
}

Vec<ZZ_pE> getExceptionalSubset(long size) {
    Vec<ZZ_pE> elms;
    elms.SetLength(size);

    for (long i = 0; i < size; i++) {
        elms[i] = indexedElementInExceptionalSet(i);
    }

    return elms;
}

Vec<ZZ_pEX> getTargetPolynomialTerms(long size) {
    Vec<ZZ_pEX> elms;
    elms.SetLength(size);

    for (long i = 0; i < size; i++) {
        std::cout << i << std::endl;
        ZZ_pEX term;
        SetX(term);
        SetCoeff(term, 0, -indexedElementInExceptionalSet(i));
        elms[i] = term;
    }

    return elms;
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

// Element must be invertible.
ZZ_pE getInverse(ZZ_pE element) {
    ZZX elemX = to_ZZX(rep(element));
    ZZX mod = to_ZZX(ZZ_pE::modulus());
    ZZX s, t;
    ZZ r;
    XGCD(r, s, t, elemX, mod, 1);

    ZZ_pX rPX = ZZ_pX();
    SetCoeff(rPX, 0, to_ZZ_p(r));
    
    ZZ_pX inverse;
    divide(inverse, to_ZZ_pX(s), rPX);

    return to_ZZ_pE(inverse); 
}