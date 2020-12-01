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

    cout << "-------Testing small example-------\n";
    {
        ZZ_pE g_1, g_2;
        g_1 = indexedElementInExceptionalSet(3);
        g_2 = indexedElementInExceptionalSet(13);
        ZZ_pE diff = g_1 - g_2;
        cout << "g_1: " << g_1 << "\n";
        cout << "g_2: " << g_2 << "\n";
        cout << "g_1-g_2: " << diff << "\n";
        cout << "g_2-g_1: " << -diff << "\n";

        ZZ_pE diffInv = getInverse(diff);
        cout << "(g1-g2)^-1: " << diffInv << "\n";
        cout << "(g1-g2) * (g1-g2)^-1: " << diff * diffInv << "\n";

        ZZ_pE minusDiffInv = getInverse(-diff);
        cout << "(g2-g1)^-1: " << minusDiffInv << "\n";
        cout << "(g2-g1) * (g2-g1)^-1: " << (-diff) * minusDiffInv << "\n";
    }
    
    cout << "-------Testing all elements in exceptional set-------\n";

    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            if (i == j) {
                continue;
            }
            ZZ_pE diff = indexedElementInExceptionalSet(i) - indexedElementInExceptionalSet(j);
            ZZ_pE diffInv = getInverse(diff);
            cout << "(" << i << "," << j << "): " << diff * diffInv << "\n";
        }
    }
}