#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ_pE.h>
#include <NTL/ZZX.h>
#include <vector>

#include <gr.h>
#include <qrp.h>
#include <setup.h>
#include <rinocchio.h>

using namespace std;
using namespace NTL;

int main() {

    ZZ modulus = ZZ(1) << 64;
    ZZ_p::init(modulus);

    // P = x^4 + x + 1
    ZZ_pX P = ZZ_pX();
    ZZ_p one = ZZ_p(1);
    SetCoeff(P, 0, one);
    SetCoeff(P, 1, one);
    SetCoeff(P, 4, one);

    // instantiate GF(2^64, 4)
    ZZ_pE::init(P);
    cout << "modulus: " << P << "\n";

    ZZ_pE g_1, g_2;
    g_1 = indexedElementInExceptionalSet(3);
    g_2 = indexedElementInExceptionalSet(13);
    ZZ_pE diff = g_1 - g_2;
    cout << "g_1: " << g_1 << "\n";
    cout << "g_2: " << g_2 << "\n";
    cout << "g_1-g_2: " << diff << "\n";
    cout << "g_2-g_1: " << -diff << "\n";

    //   XGCD(d, x, t, a, f); // x result, a elem, f mod, d, t fresh
    // ZZ_pX dx, tx, x;
    // ZZ_pX dixx = rep(diff);
    // PlainXGCD(dx, x, tx, dixx, ZZ_pE::modulus());
    // cout << x << "inv\n";

    // cout << "lead coeff g_1-g_2: " << LeadCoeff(rep(diff)) << "\n";
    // cout << "lead coeff g_2-g_1: " << LeadCoeff(rep(-diff)) << "\n";
    
    {
        ZZX diffx = to_ZZX(rep(diff));
        ZZX q, r;
        // HomPseudoDivRem(q, r, to_ZZX(ZZ_pE::modulus()), diffx);
        ZZX inverse, mod, t;
        mod = to_ZZX(ZZ_pE::modulus());
        ZZ d;
        XGCD(d, inverse, t, diffx, mod, 1);
        cout << "d % 2^64: " << d % modulus << "\n";
        ZZ_p dP = to_ZZ_p(d);
        ZZ_pX dPX = ZZ_pX();
        SetCoeff(dPX, 0, dP);

        ZZ_pX invX = to_ZZ_pX(inverse);
        divide(invX, invX, dPX);
        ZZ_pE actualInverse = to_ZZ_pE(invX);
        cout << "Inverse: " << actualInverse << "\n";

        ZZ_pE reducedOne = actualInverse * diff;
        cout << "inverse * diff: " <<  reducedOne  << "\n"; 
        cout << "g_1-g_2 passed\n";
    }
}