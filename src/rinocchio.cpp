#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>


using namespace std;
using namespace NTL;

ZZ_pE randomNonZeroInExceptionalSet() {//TODO
    // d random bit coefficients
    // use IsZero()
    return conv<ZZ_pE>("[1 1 1 1]"); // 1 + x + x^2 + x^3?
}

ZZ_pE randomInvertible() {
    while (true) {
        // random element in R:
        ZZ_pE res = random_ZZ_pE();
        cout << rep(res) << "\n";
        //When d sufficiently large 
        //Check at least one coefficient is one
        for (int i = 0; i < ZZ_pE::degree(); i++) {
            if (!divide(rep(rep(res)[i]), 2)) return res;
        }
    }
}

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

    //Find non-zero s in exceptional set (i.e. in A^*):
    ZZ_pE s = randomNonZeroInExceptionalSet();
    cout << rep(s) << "\n" ;
    cout << s << "\n" ;
    //Find random r_v, r_w in R^*, i.e. d random coefficients where at least one is odd
    ZZ_pE r_v, r_w, r_y;
    r_v = randomInvertible();
    r_w = randomInvertible();
    r_y = r_v * r_w;
    ZZ_pE alpha, alpha_v, alpha_w, alpha_y;
    alpha = randomInvertible();
    alpha_v = randomInvertible();
    alpha_w = randomInvertible();
    alpha_y = randomInvertible();

    ZZ_pE beta;
    do {
        beta = random_ZZ_pE();
    } while (IsZero(beta));

    //TODO: keypair
    //TODO: QRP
    //TODO: CRS
}