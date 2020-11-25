#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ_pE.h>
#include <vector>


using namespace std;
using namespace NTL;

//dummy encryption:
ZZ_pE E(ZZ_pE x) {return x;}

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

// Random element in A
ZZ_pE randomInExceptionalSet() {
    ZZ_pX a = ZZ_pX();
    for (int i = 0; i < ZZ_pE::degree(); i++) {
        long coeff;
        RandomBits(coeff, 1);
        SetCoeff(a, i, coeff);
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

    //TODO: keypair, when joye-libert is working
    
    
    //QRP:

    //Pick distinct elements of exceptional set for each gate:
    ZZ_pE g_1, g_2;
    g_1 = randomInExceptionalSet();
    do {
        g_2 = randomInExceptionalSet();
    } while (g_1 == g_2);

    // Compute t(x) = (x - g_1) * (x - g_2)
    Vec<ZZ_pEX> V, W, Y;
    V.SetLength(6);
    W.SetLength(6);
    Y.SetLength(6);
    
    ZZ_pEX t;
    {
        ZZ_pEX x;
        SetX(x);
        t = (x - g_1)*(x-g_2);
        cout << "g_1" << g_1 << "\n";
        cout << "g_2" << g_2 << "\n";
        cout << "x" << x << "\n";
        cout << "t" << t << "\n";    
        cout << "t(g_1)" << eval(t, g_1) << "\n";    
        cout << "t(g_2)" << eval(t, g_2) << "\n";    
        cout << "t(g_2+g_1)" << eval(t, g_2+g_1) << "\n";    

        ZZ_pE galloisOne = conv<ZZ_pE>(1);

        Vec<ZZ_pE> a;
        a.append(g_1);
        a.append(g_2);
        Vec<ZZ_pE> b_v_1;
        b_v_1.append(ZZ_pE::zero());
        b_v_1.append(galloisOne);
        Vec<ZZ_pE> b_v_2;
        b_v_2.append(ZZ_pE::zero());
        b_v_2.append(galloisOne);
        Vec<ZZ_pE> b_v_3;
        b_v_3.append(galloisOne);
        b_v_3.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_v_4;
        b_v_4.append(ZZ_pE::zero());
        b_v_4.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_v_5;
        b_v_5.append(ZZ_pE::zero());
        b_v_5.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_v_6;
        b_v_6.append(ZZ_pE::zero());
        b_v_6.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_w_1;
        b_w_1.append(ZZ_pE::zero());
        b_w_1.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_w_2;
        b_w_2.append(ZZ_pE::zero());
        b_w_2.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_w_3;
        b_w_3.append(ZZ_pE::zero());
        b_w_3.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_w_4;
        b_w_4.append(galloisOne);
        b_w_4.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_w_5;
        b_w_5.append(ZZ_pE::zero());
        b_w_5.append(galloisOne);
        Vec<ZZ_pE> b_w_6;
        b_w_6.append(ZZ_pE::zero());
        b_w_6.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_1;
        b_y_1.append(ZZ_pE::zero());
        b_y_1.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_2;
        b_y_2.append(ZZ_pE::zero());
        b_y_2.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_3;
        b_y_3.append(ZZ_pE::zero());
        b_y_3.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_4;
        b_y_4.append(ZZ_pE::zero());
        b_y_4.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_5;
        b_y_5.append(galloisOne);
        b_y_5.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_6;
        b_y_6.append(ZZ_pE::zero());
        b_y_6.append(galloisOne);


        interpolate(V(1), a, b_v_1);
        interpolate(V(2), a, b_v_2);
        interpolate(V(3), a, b_v_3);
        interpolate(V(4), a, b_v_4);
        interpolate(V(4), a, b_v_4);
        interpolate(V(5), a, b_v_5);
        interpolate(W(6), a, b_w_6);
        interpolate(W(1), a, b_w_1);
        interpolate(W(2), a, b_w_2);
        interpolate(W(3), a, b_w_3);
        interpolate(W(4), a, b_w_4);
        interpolate(W(5), a, b_w_5);
        interpolate(W(6), a, b_w_6);
        interpolate(Y(1), a, b_y_1);
        interpolate(Y(2), a, b_y_2);
        interpolate(Y(3), a, b_y_3);
        interpolate(Y(4), a, b_y_4);
        interpolate(Y(5), a, b_y_5);
        interpolate(Y(6), a, b_y_6);
    }
    
    //TODO: CRS
    //{E(S^i)}_i=0^d
    //{E(alpha * S^i)}_i=0^d
    Vec<ZZ_pE> powersOfS;
    Vec<ZZ_pE> powersOfSMultAlpa;
    {
        ZZ_pE acc;
        set(acc);
        for (int i = 0; i <= ZZ_pE::degree(); i++) {
            ZZ_pE ithPower = ZZ_pE(acc);
            powersOfS.append(ithPower);

            ZZ_pE ithPowerMultAlpha;
            mul(ithPowerMultAlpha, alpha, ithPower);
            powersOfSMultAlpa.append(ithPowerMultAlpha);
            
            mul(acc, acc, s);
        }
    }

    
    Vec<ZZ_pE> rvVofS, rwWofS, ryYofS, alpharvVofS, alpharwWofS, alpharyYofS, betaSums;
    rvVofS.SetLength(6);
    rwWofS.SetLength(6);
    ryYofS.SetLength(6);
    alpharvVofS.SetLength(6);
    alpharwWofS.SetLength(6);
    alpharyYofS.SetLength(6);
    betaSums.SetLength(6);

    for (int k = 0; k < V.length(); k++) {
        //{E(r_v * v_k(S))}_k\in I_mid'
        //{E(r_w * w_k(S))}_k\in I_mid
        //{E(r_y * y_k(S))}_k\in I_mid
        ZZ_pE rvVkofS, rwWkofS, ryYkofS;
        eval(rvVkofS, V[k], s);
        eval(rwWkofS, W[k], s);
        eval(ryYkofS, Y[k], s);
        mul(rvVkofS, rvVkofS, r_v);
        mul(rwWkofS, rwWkofS, r_w);
        mul(ryYkofS, ryYkofS, r_y);
        rvVofS[k] = E(rvVkofS);
        rwWofS[k] = E(rwWkofS);
        ryYofS[k] = E(ryYkofS);
            
        //{E(alpha_v * r_v * v_k(S))}_k\in I_mid
        //{E(alpha_w * r_w * w_k(S))}_k\in I_mid
        //{E(alpha_y * r_y * y_k(S))}_k\in I_mid
        ZZ_pE alpharvVkofS, alpharwWkofS, alpharyYkofS;
        mul(alpharvVkofS, rvVkofS, alpha_v);
        mul(alpharwWkofS, rwWkofS, alpha_w);
        mul(alpharyYkofS, ryYkofS, alpha_y);
        alpharvVofS[k] = E(alpharvVkofS);
        alpharwWofS[k] = E(alpharwWkofS);
        alpharyYofS[k] = E(alpharyYkofS);

        //{E(beta ((r_v * v_k(S)) + (r_w * w_k(S)) + (r_y * y_k(S)))}_k\in I_mid
        ZZ_pE kthBetaSum;
        kthBetaSum = beta * (alpharvVkofS + alpharwWkofS + alpharyYkofS);
        betaSums[k] = E(kthBetaSum);
    }

    //pk
}