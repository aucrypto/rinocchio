#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ_pE.h>
#include <NTL/ZZX.h>
#include <vector>
#include <NTL/vec_vec_ZZ_p.h>

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
    cout << "g_1: " << g_1 << "\n";
    Vec<ZZ_pE> multWires;
    multWires.append(g_1);
    multWires.append(g_2);

    Vec<ZZ_pE> v_values;
    ZZ_pE galloisOne;
    set(galloisOne);
    v_values.append(galloisOne);
    v_values.append(ZZ_pE::zero());

    ZZ_pEX poly;
    cout << "CHECKPOINT 1\n";
    
    ZZ_pEX delta_g1, delta_g2;
    ZZ_pEX x = ZZ_pEX();
    SetX(x);

    set(delta_g1);
    set(delta_g2);
    cout << delta_g2 << "\n";

    delta_g1 = x - g_2;
    mul(delta_g1, delta_g1, getInverse(g_1 - g_2));

    delta_g2 = x - g_1;
    cout << delta_g2 << "\n";
    mul(delta_g2, delta_g2, getInverse(g_2 - g_1));
    cout << delta_g2 << "\n";

    poly = delta_g2; 
    
    cout << "poly(g1): " << eval(poly, g_1) << "\n";
    cout << "poly(g2): " << eval(poly, g_2) << "\n";

    return 0;
}


