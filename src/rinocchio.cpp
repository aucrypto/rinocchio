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

CRS getCRS(QRP prog, secretState ss) {
    
    //{E(S^i)}_i=0^d
    //{E(alpha * S^i)}_i=0^d
    Vec<Vec<ZZ>> powersOfS;
    Vec<Vec<ZZ>> powersOfSMultAlpha;
    {
        ZZ_pE acc;
        set(acc);
        for (int i = 0; i <= prog.circuit.numberOfMultiplicationGates; i++) {//todo == d?
            ZZ_pE ithPower = ZZ_pE(acc);
            powersOfS.append(E(ithPower));

            ZZ_pE ithPowerMultAlpha;
            mul(ithPowerMultAlpha, ss.alpha, ithPower);
            powersOfSMultAlpha.append(E(ithPowerMultAlpha));
            
            mul(acc, acc, ss.s);
        }
    }

    long sizeOfImid = prog.circuit.numberOfMidWires;
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
        eval(rvVkofS, prog.V[k+prog.midOffset], ss.s);
        eval(rwWkofS, prog.W[k+prog.midOffset], ss.s);
        eval(ryYkofS, prog.Y[k+prog.midOffset], ss.s);
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

Proof prove(QRP prog, CRS crs, Vec<ZZ_p> allWireValues) {


    // Compute p = V*W-Y
    // P = W * W * Y = (Sum c_k * v_k(x)) * (Sum c_k * w_k(x)) - (Sum c_k * y_k(x))
    ZZ_pEX V, W, Y;
    for (int k = 0; k < prog.circuit.numberOfWires; k++) {
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
    Vec<ZZ> rvVmidOfS;
    Vec<ZZ> rwWmidOfS;
    Vec<ZZ> ryYmidOfS;
    Vec<ZZ> alphaVrvVmidOfS;
    Vec<ZZ> alphaWrwWmidOfS;
    Vec<ZZ> alphaYryYmidOfS;
    Vec<ZZ> betaSum;
    rvVmidOfS.SetLength(ZZ_pE::degree());
    rwWmidOfS.SetLength(ZZ_pE::degree());
    ryYmidOfS.SetLength(ZZ_pE::degree());
    alphaVrvVmidOfS.SetLength(ZZ_pE::degree());
    alphaWrwWmidOfS.SetLength(ZZ_pE::degree());
    alphaYryYmidOfS.SetLength(ZZ_pE::degree());
    betaSum.SetLength(ZZ_pE::degree());
    for (int k = prog.midOffset; k < prog.outOffset; k++) {
        rvVmidOfS += rep(allWireValues[k]) * crs.rvVofS[k - prog.midOffset];
        rwWmidOfS += rep(allWireValues[k]) * crs.rwWofS[k - prog.midOffset];
        ryYmidOfS += rep(allWireValues[k]) * crs.ryYofS[k - prog.midOffset];
        alphaVrvVmidOfS += rep(allWireValues[k]) * crs.alpharvVofS[k - prog.midOffset];
        alphaWrwWmidOfS += rep(allWireValues[k]) * crs.alpharwWofS[k - prog.midOffset];
        alphaYryYmidOfS += rep(allWireValues[k]) * crs.alpharyYofS[k - prog.midOffset];
        betaSum += rep(allWireValues[k]) * crs.betaSums[k - prog.midOffset];

    }

    // E(h(s))
    // E(alpha * h(s))
    ZZX hOfS, alphaHofS;
    Vec<ZZ> mod = to_vec_ZZ(ZZ_pE::modulus().f.rep);
    ZZX modX = to_ZZX(mod);
    for (int i = 0 ; i <= deg(H); i++) {
        //todo How do we know that the degree of H is at most the number of multiplication gates?
        Vec<ZZ> ithCoeffOfH = to_vec_ZZ(rep(coeff(H, i)).rep);
        ZZX ithCoeffOfHX, ithPowerOfSX, alphaIthPowerOfSX;
        ithCoeffOfHX = to_ZZX(ithCoeffOfH);
        ithPowerOfSX = to_ZZX(crs.powersOfS[i]);
        alphaIthPowerOfSX = to_ZZX(crs.powersOfSMultAlpha[i]);
        ZZX ithTermOfHSX = (ithCoeffOfHX * ithPowerOfSX) % modX;
        ZZX ithTermOfalphaHSX = (ithCoeffOfHX * alphaIthPowerOfSX) % modX;
        hOfS += ithTermOfHSX;
        alphaHofS += ithTermOfalphaHSX;
    }
    hOfS %= modX;//todo do more reductions
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
    if (secret.alpha * hOfS != alphaHOfS) {
        cout << "\n";
        return false;
    }

    // todo compute P_io
    // compute P_in
    ZZ_pE v_io, w_io, y_io;
    for (int k = 0; k < input.length(); k++) {
        v_io += input[k] * eval(qrp.V[k], secret.s);
        w_io += input[k] * eval(qrp.W[k], secret.s);
        y_io += input[k] * eval(qrp.Y[k], secret.s);
    }
    // compute P_out
    for (int k = 0; k < output.length(); k++) {
        v_io += output[k] * eval(qrp.V[k+qrp.outOffset], secret.s);
        w_io += output[k] * eval(qrp.W[k+qrp.outOffset], secret.s);
        y_io += output[k] * eval(qrp.Y[k+qrp.outOffset], secret.s);
    }

    //todo we cannot divide by r_v, r_w, and r_y because of the inv bug, but we have r_y = r_v * r_w
    ZZ_pE computedV = secret.r_v * v_io + rvVmidOfS;
    ZZ_pE computedW = secret.r_w * w_io + rwWmidOfS;
    ZZ_pE computedY = secret.r_y * y_io + ryYmidOfS;

    ZZ_pE computedP = computedV * computedW - computedY;
    ZZ_pE hMultT = hOfS * eval(qrp.t, secret.s);
    if (computedP != secret.r_y * hMultT) {
        cout << computedP << "\n";
        return false;
    }

    return  true;
}
