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

using namespace std;
using namespace NTL;

//dummy encryption:
Vec<ZZ> E(ZZ_pE x) {
    Vec<ZZ> res;
    // res = to_vec_ZZ(rep(x).rep);
    for (int i = 0; i < ZZ_pE::degree(); i++) {
        ZZ_p coeff;
        GetCoeff(coeff, rep(x), i);//returns 0 if i out of range.
        //Todo actually encrypt the coefficient under the ith pulic key.
        res.append(rep(coeff));
    }
    return res;
}

// dummy decryption
ZZ_pE D(Vec<ZZ> y) {
    ZZ_pX x;
    for (int i = 0; i < y.length(); i++) {
        SetCoeff(x, i, conv<ZZ_p>(y[i]));
    }
    return conv<ZZ_pE>(x);
}

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
        kthBetaSum = ss.beta * (rvVkofS + rwWkofS + ryYkofS);
        betaSums[k] = E(kthBetaSum);
    }
    //pk

    CRS crs;
    crs.rvVofS = rvVofS;
    crs.rwWofS = rwWofS;
    crs.ryYofS = ryYofS;
    crs.alpharvVofS = alpharvVofS;
    crs.alpharwWofS = alpharwWofS;
    crs.alpharyYofS = alpharyYofS;
    crs.betaSums = betaSums;
    crs.powersOfS = powersOfS;
    crs.powersOfSMultAlpha = powersOfSMultAlpha;
    return crs;
}

struct Proof {
    Vec<ZZ> rvVmidOfS;
    Vec<ZZ> rwWmidOfS;
    Vec<ZZ> ryYmidOfS;
    Vec<ZZ> alphaVrvVmidOfS;
    Vec<ZZ> alphaWrwWmidOfS;
    Vec<ZZ> alphaYryYmidOfS;
    Vec<ZZ> betaSum;
    Vec<ZZ> hOfS;
    Vec<ZZ> alphaHOfS;
};

Proof prove(QRP prog, CRS crs, Vec<ZZ_p> input) {
    Vec<ZZ_p> allWireValues = input;
    //Compute values for all multiplication gates:
    ZZ_p c_5, c_6;//todo generalize using circuit representation
    
    c_5 = input(3) * input(4);
    allWireValues.append(c_5);
    cout << c_5 << " should be 6\n";

    c_6 = c_5 * (input(1) + input(2));
    allWireValues.append(c_6);
    cout << c_6 << " should be 42\n";

    // Compute p = V*W-Y
    // P = W * W * Y = (Sum c_k * v_k(x)) * (Sum c_k * w_k(x)) - (Sum c_k * y_k(x))
    ZZ_pEX V, W, Y;
    for (int k = 0; k < prog.numberOfWires; k++) {
        V += allWireValues[k] * prog.V[k];
        W += allWireValues[k] * prog.W[k];
        Y += allWireValues[k] * prog.Y[k];
    }
    ZZ_pEX P = V*W-Y;

    // Compute h = p*(t^-1)
    ZZ_pEX H = P / prog.t;//todo we should precompute the inverse of t instead.


    // E(r_v * Vmid(S))
    // E(r_v * Vmid(S) * alpha_v)
    // E(r_w * Wmid(S))
    // E(r_w * Wmid(S) * alpha_w)
    // E(r_y * Ymid(S))
    // E(r_y * Ymid(S) * alpha_y)
    // E(beta( (r_v * Vmid(S)) + (r_w * Wmid(S)) +(r_y * Ymid(S)) ))

    // for now the Vec's in the crs only contain one element:
    Vec<ZZ> rvVmidOfS = rep(c_5) * crs.rvVofS[0];
    Vec<ZZ> rwWmidOfS = rep(c_5) * crs.rwWofS[0];
    Vec<ZZ> ryYmidOfS = rep(c_5) * crs.ryYofS[0];
    Vec<ZZ> alphaVrvVmidOfS = rep(c_5) * crs.alpharvVofS[0];
    Vec<ZZ> alphaWrwWmidOfS = rep(c_5) * crs.alpharwWofS[0];
    Vec<ZZ> alphaYryYmidOfS = rep(c_5) * crs.alpharyYofS[0];
    Vec<ZZ> betaSum = rep(c_5) * crs.betaSums[0];

    // E(h(s))
    // E(alpha * h(s))

    // todo: `evaluate' h(S) by mutiplying powers of S by coefficients of h
    // todo: `evaluate' alpha*h(S) by mutiplying powers of alpha*S by coefficients of h
    ZZX hOfS, alphaHofS;
    Vec<ZZ> mod = to_vec_ZZ(ZZ_pE::modulus().f.rep);
    ZZX modX = to_ZZX(mod);
    for (int i = 0 ; i < ZZ_pE::degree(); i++) {//todo d + 1?
        Vec<ZZ> ithCoeffOfH = to_vec_ZZ(rep(coeff(H, i)).rep);// todo don't use the actual encryption
        // I guess we need to multiply them as polynomials? Use ZZX?
        ZZX ithCoeffOfHX, ithPowerOfSX, alphaIthPowerOfSX; //todo only makes sense if ciphertexts are "straightforwardly" multiplied together..
        ithCoeffOfHX = to_ZZX(ithCoeffOfH);
        ithPowerOfSX = to_ZZX(crs.powersOfS[i]);
        alphaIthPowerOfSX = to_ZZX(crs.powersOfSMultAlpha[i]);
        ZZX ithTermOfHSX = (ithCoeffOfHX * ithPowerOfSX) % modX;
        ZZX ithTermOfalphaHSX = (ithCoeffOfHX * alphaIthPowerOfSX) % modX;
        hOfS += ithTermOfHSX;
        alphaHofS += ithTermOfalphaHSX;//todo check
    }
    hOfS %= modX;
    alphaHofS %= modX;

    return Proof{
        .rvVmidOfS = rvVmidOfS,
        .rwWmidOfS = rwWmidOfS,
        .ryYmidOfS  =ryYmidOfS,
        .alphaVrvVmidOfS = alphaVrvVmidOfS,
        .alphaWrwWmidOfS = alphaWrwWmidOfS,
        .alphaYryYmidOfS = alphaYryYmidOfS,
        .betaSum = betaSum,
        .hOfS = hOfS.rep,
        .alphaHOfS = alphaHofS.rep,
    };
}

bool verify(QRP qrp, secretState secret, CRS crs, Proof pi, Vec<ZZ_p> input, Vec<ZZ_p> output) {
    ZZ_pE rvVmidOfS = D(pi.rvVmidOfS);
    ZZ_pE rwWmidOfS = D(pi.rwWmidOfS);
    ZZ_pE ryYmidOfS = D(pi.ryYmidOfS);
    ZZ_pE alphavrvVmidOfS = D(pi.alphaVrvVmidOfS);
    ZZ_pE alphawrwWmidOfS = D(pi.alphaWrwWmidOfS);
    ZZ_pE alphayryYmidOfS = D(pi.alphaYryYmidOfS);
    if (secret.alpha_v * rvVmidOfS != alphavrvVmidOfS) {
        cout << "alpha_v * r_v * v_mid(S) != \n";
        return false;
    }
    if (secret.alpha_w * rwWmidOfS != alphawrwWmidOfS) {
        cout << "alpha_w * r_w * w_mid(S) != \n";
        return false;
    }
    if (secret.alpha_y * ryYmidOfS != alphayryYmidOfS) {
        cout << "alpha_y * r_y * y_mid(S) != \n";
        return false;
    }
    ZZ_pE betaSum = D(pi.betaSum);
    ZZ_pE computedBetaSum = secret.beta * (rvVmidOfS + rwWmidOfS + ryYmidOfS);

    if (computedBetaSum != betaSum) {
        cout << secret.beta << "L != \n";
        return false;
    }
    
    ZZ_pE hOfS = D(pi.hOfS);
    ZZ_pE alphaHOfS = D(pi.alphaHOfS);



    return  true;
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

    Vec<ZZ_p> output;
    output.append(ZZ_p(42));

    cout << verify(qrp, state, crs, pi, input, output) << "\n";
}