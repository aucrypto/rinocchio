#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>


using namespace std;
using namespace NTL;

ZZ_pE randomInvertible() {
    while (true) {
        // random element in R:
        ZZ_pE res = random_ZZ_pE();
        //When d sufficiently large 
        //Check at least one coefficient is one
        for (int i = 0; i < ZZ_pE::degree(); i++) {
            ZZ coeff = rep(rep(res)[i]);//todo remove tmp?
            if (IsOdd(coeff)) return res;
            
        }
    }
}

ZZ_pE randomNonZeroInExceptionalSet() {//TODO it's not pretty but it works for now :D
    ZZ_pE preMod = randomInvertible(); //i.e. not zero when coefficients reduced mod 2
    // bool reducedCoeffs[ZZ_pE::degree()]; //todo, why does this work? Doesn't it need to be a constant?
    ZZ_pX reducedModTwo = ZZ_pX();
    for (int i = 0; i < ZZ_pE::degree(); i++) {
        // bool reducedCoeff = !divide(rep(rep(preMod)[i]), 2);
        // reducedCoeffs[i] = reducedCoeff;
        SetCoeff(reducedModTwo, i, rem(rep(rep(preMod)[i]), 2));
    }
    // ZZ_pE fromarray = conv<ZZ_pE>(reducedCoeffs);
    ZZ_pE fromPX = to_ZZ_pE(reducedModTwo);
    return fromPX;
    // d random bit coefficients
    // use IsZero()
    // return conv<ZZ_pE>("[1 1 1 1]"); // 1 + x + x^2 + x^3?
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