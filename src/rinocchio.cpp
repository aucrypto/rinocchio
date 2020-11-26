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
Vec<ZZ> E(ZZ_pE x) {
    Vec<ZZ> res;
    ZZ_pX polynomialRep = rep(x);
    for (int i = 0; i < ZZ_pE::degree(); i++) {
        ZZ_p coeff;
        GetCoeff(coeff, polynomialRep, i);//returns 0 if i out of range.
        //Todo actually encrypt the coefficient under the ith pulic key.
        res.append(rep(coeff));
    }
    return res;
}

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

struct QRP {
    long numberOfWires;
    long numberOfInputWires;
    long numberOfOutputWires;
    ZZ_pEX t; //target polynomial
    Vec<ZZ_pEX> V;
    Vec<ZZ_pEX> W;
    Vec<ZZ_pEX> Y;
};

struct CRS {
    Vec<Vec<ZZ>> powersOfS;
    Vec<Vec<ZZ>> powersOfSMultAlpha;
    Vec<Vec<ZZ>> rvVofS;
    Vec<Vec<ZZ>> rwWofS;
    Vec<Vec<ZZ>> ryYofS;
    Vec<Vec<ZZ>> alpharvVofS;
    Vec<Vec<ZZ>> alpharwWofS;
    Vec<Vec<ZZ>> alpharyYofS;
    Vec<Vec<ZZ>> betaSums; 
    bool publicKey;
};

QRP getQRP() {
    //Pick distinct elements of exceptional set for each gate:
    ZZ_pE g_1, g_2;
    g_1 = indexedElementInExceptionalSet(1);
    g_2 = indexedElementInExceptionalSet(2);
    ZZ_pE diff = g_1 - g_2;
    // cout << "g_1" << g_1 << "\n";
    // cout << "g_2" << g_2 << "\n";
    // cout << "g_1-g_2" << diff << "\n";
    // cout << "g_2-g_1" << -diff << "\n";
    // // cout << "lead coeff g_1-g_2: " << LeadCoeff(rep(diff)) << "\n";
    // // cout << "lead coeff g_2-g_1: " << LeadCoeff(rep(-diff)) << "\n";
    // inv(diff);
    // cout << "g_1-g_2 passed\n";
    // inv(-diff);
    // cout << "g_2-g_1 passed\n";



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
        // cout << "g_1" << g_1 << "\n";
        // cout << "g_2" << g_2 << "\n";
        // cout << "x" << x << "\n";
        // cout << "t" << t << "\n";    
        // cout << "t(g_1)" << eval(t, g_1) << "\n";    
        // cout << "t(g_2)" << eval(t, g_2) << "\n";    
        // cout << "t(g_2+g_1)" << eval(t, g_2+g_1) << "\n";    

        // ZZ_pE galloisOne = conv<ZZ_pE>(1);
        ZZ_pE galloisOne;
        set(galloisOne);

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
        interpolate(V(5), a, b_v_5);
        interpolate(V(6), a, b_v_6);

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
    
    QRP qrp;
    qrp.numberOfWires = 6;
    qrp.numberOfInputWires = 4;
    qrp.numberOfOutputWires = 1;
    
    qrp.t  = t;
    qrp.V = V;
    qrp.W = W;
    qrp.Y = Y;
    return qrp;
}


struct secretState {
    ZZ_pE s, r_v, r_w, r_y, alpha, alpha_v, alpha_w, alpha_y, beta;
    bool secretKey;
};

secretState setup() {
    //Find non-zero s in exceptional set (i.e. in A^*):
    ZZ_pE s;
    s = randomNonZeroInExceptionalSet();
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

    secretState ss;
    ss.s = s;
    ss.r_v = r_v;
    ss.r_w = r_w;
    ss.r_y = r_y;
    ss.alpha = alpha;
    ss.alpha_v = alpha_v;
    ss.alpha_w = alpha_w;
    ss.alpha_y = alpha_y;
    ss.beta = beta;
    //TODO: keypair, when joye-libert is working
    
    ss.secretKey = false;
    
    return ss;
}

CRS getCRS(QRP prog, secretState ss) {
    
    //{E(S^i)}_i=0^d
    //{E(alpha * S^i)}_i=0^d
    Vec<Vec<ZZ>> powersOfS;
    Vec<Vec<ZZ>> powersOfSMultAlpha;
    {
        ZZ_pE acc;
        set(acc);
        for (int i = 0; i < ZZ_pE::degree(); i++) {//todo == d?
            ZZ_pE ithPower = ZZ_pE(acc);
            powersOfS.append(E(ithPower));

            ZZ_pE ithPowerMultAlpha;
            mul(ithPowerMultAlpha, ss.alpha, ithPower);
            powersOfSMultAlpha.append(E(ithPowerMultAlpha));
            
            mul(acc, acc, ss.s);
        }
    }

    long sizeOfImid = prog.numberOfWires - prog.numberOfInputWires - prog.numberOfOutputWires;
    Vec<Vec<ZZ>> rvVofS, rwWofS, ryYofS, alpharvVofS, alpharwWofS, alpharyYofS, betaSums;
    rvVofS.SetLength(sizeOfImid);
    rwWofS.SetLength(sizeOfImid);
    ryYofS.SetLength(sizeOfImid);
    alpharvVofS.SetLength(sizeOfImid);
    alpharwWofS.SetLength(sizeOfImid);
    alpharyYofS.SetLength(sizeOfImid);
    betaSums.SetLength(sizeOfImid);

    for (int k = 0; k < sizeOfImid; k++) {
        //{E(r_v * v_k(S))}_k\in I_mid'
        //{E(r_w * w_k(S))}_k\in I_mid
        //{E(r_y * y_k(S))}_k\in I_mid
        ZZ_pE rvVkofS, rwWkofS, ryYkofS;
        eval(rvVkofS, prog.V[k+prog.numberOfInputWires], ss.s);
        eval(rwWkofS, prog.W[k+prog.numberOfInputWires], ss.s);
        eval(ryYkofS, prog.Y[k+prog.numberOfInputWires], ss.s);
        mul(rvVkofS, rvVkofS, ss.r_v);
        mul(rwWkofS, rwWkofS, ss.r_w);
        mul(ryYkofS, ryYkofS, ss.r_y);
        rvVofS[k] = E(rvVkofS);
        rwWofS[k] = E(rwWkofS);
        ryYofS[k] = E(ryYkofS);
            
        //{E(alpha_v * r_v * v_k(S))}_k\in I_mid
        //{E(alpha_w * r_w * w_k(S))}_k\in I_mid
        //{E(alpha_y * r_y * y_k(S))}_k\in I_mid
        ZZ_pE alpharvVkofS, alpharwWkofS, alpharyYkofS;
        mul(alpharvVkofS, rvVkofS, ss.alpha_v);
        mul(alpharwWkofS, rwWkofS, ss.alpha_w);
        mul(alpharyYkofS, ryYkofS, ss.alpha_y);
        alpharvVofS[k] = E(alpharvVkofS);
        alpharwWofS[k] = E(alpharwWkofS);
        alpharyYofS[k] = E(alpharyYkofS);

        //{E(beta ((r_v * v_k(S)) + (r_w * w_k(S)) + (r_y * y_k(S)))}_k\in I_mid
        ZZ_pE kthBetaSum;
        kthBetaSum = ss.beta * (alpharvVkofS + alpharwWkofS + alpharyYkofS);
        betaSums[k] = E(kthBetaSum);
    }
    //pk

    CRS crs;
    return crs;
}

struct Proof {
    bool isTrue;

};

Proof prove(QRP prog, CRS crs, Vec<ZZ_p> input) {
    //Compute values for all multiplication gates:
    ZZ_p c_5, c_6;
    c_5 = input(3) * input(4);
    cout << c_5 << " should be 6\n";
    c_6 = c_5 * (input(1) + input(2));
    cout << c_6 << " should be 42\n";

    // E(r_v * Vmid(S))
    // E(r_v * Vmid(S) * alpha_v)
    // E(r_w * Wmid(S))
    // E(r_w * Wmid(S) * alpha_w)
    // E(r_y * Ymid(S))
    // E(r_y * Ymid(S) * alpha_y)
    // E(beta( (r_v * Vmid(S)) + (r_w * Wmid(S)) +(r_y * Ymid(S)) ))

    // ZZ_p rvVmidOfS = c_5 * crs.rvVofS;
    ZZ_p alphavrvVmidOfS;
    ZZ_p rwWmidOfS;
    ZZ_p alphawrwWmidOfS;
    ZZ_p ryYmidOfS;
    ZZ_p alphayryYmidOfS;
    ZZ_p betaSum;

    // E(h(s))
    // E(alpha * h(s))

    Proof proof;
    return proof;
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

    QRP qrp = getQRP();
    secretState state = setup();
    CRS crs = getCRS(qrp, state);
    Vec<ZZ_p> input;
    input.append(ZZ_p(3));
    input.append(ZZ_p(4));
    input.append(ZZ_p(2));
    input.append(ZZ_p(3));
    Proof pi = prove(qrp, crs, input);
}