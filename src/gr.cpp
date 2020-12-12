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


ZZ_pX primitiveIrredPoly(long degree) {
    ZZ_pX P;
    if (degree < 2) {
        std::cout << "Degree of quotion polynomial was " << degree << " but must be at least 2\n";
        return P;
    }
    P = 1;
    SetCoeff(P, degree);

    switch (degree)
    {
    case 2:
        // x^2 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 3:
        // x^3 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 4:
        // x^4 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 5:
        // x^5 + x^2 + 1
        SetCoeff(P, 2);
        break;
    case 6:
        // x^6 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 7:
        // x^7 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 8:
        // x^8 + x^4 + x^3 + x^2 + 1
        SetCoeff(P, 4);
        SetCoeff(P, 3);
        SetCoeff(P, 2);
        break;
    case 9:
        // x^9 + x^4 + 1
        SetCoeff(P, 4);
        break;
    case 10:
        // x^10 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 11:
        // x^11 + x^2 + 1
        SetCoeff(P, 2);
        break;
    case 12:
        // x^12 + x^6 + x^4 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 4);
        SetCoeff(P, 6);
        break;
    case 13:
        // x^13 + x^4 + x^3 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 3);
        SetCoeff(P, 4);
        break;
    case 14:
        // x^14 + x^8 + x^6 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 6);
        SetCoeff(P, 8);
        break;
    case 15:
        // x^15 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 16:
        // x^16 + x^9 + x^8 + x^7 + x^6 + x^4 + x^3 + x^2 + 1
        SetCoeff(P, 2);
        SetCoeff(P, 3);
        SetCoeff(P, 4);
        SetCoeff(P, 6);
        SetCoeff(P, 7);
        SetCoeff(P, 8);
        SetCoeff(P, 9);
        break;
    case 17:
        // x^17 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 18:
        // x^18 + x^5 + x^4 + x^3 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 3);
        SetCoeff(P, 4);
        SetCoeff(P, 5);
        break;
    case 19:
        // x^19 + x^5 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 5);
        break;
    case 20:
        // x^20 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 21:
        // x^21 + x^2 + 1
        SetCoeff(P, 2);
        break;
    case 22:
        // x^22 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 23:
        // x^23 + x^5 + 1
        SetCoeff(P, 5);
        break;
    case 24:
        // x^24 + x^7 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 7);
        break;
    case 25:
        // x^25 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 26:
        // x^26 + x^6 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 6);
        break;
    case 27:
        // x^27 + x^5 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 5);
        break;
    case 28:
        // x^28 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 29:
        // x^29 + x^2 + 1
        SetCoeff(P, 2);
        break;
    case 30:
        // x^30 + x^23 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 23);
        break;
    case 31:
        // x^31 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 32:
        // x^32 + x^22 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 22);
        break;
    case 40: //For soundness error 2^40
        // x^40 + x^29 + x^27 + x^23 + 1
        SetCoeff(P, 23);
        SetCoeff(P, 27);
        SetCoeff(P, 29);
        break;
    case 60: //For soundness error 2^60
        // x^60 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 80: //For soundness error 2^80
        // x^80 + x^75 + x^27 + x^17 + 1
        SetCoeff(P, 17);
        SetCoeff(P, 27);
        SetCoeff(P, 75);
        break;
    default:
        std::cout << "Unsupported extension degree: " << degree << std::endl;
        break;
    }
    return P;
}
