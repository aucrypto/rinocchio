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
#include <joye_libert.h>

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
    Vec<JLEncoding> powersOfS;
    Vec<JLEncoding> powersOfSMultAlpha;
    {
        ZZ_pE acc;
        set(acc);
        for (int i = 0; i <= prog.circuit.numberOfMultiplicationGates; i++) {//todo == d?
            ZZ_pE ithPower = ZZ_pE(acc);
            powersOfS.append(encode(ithPower, ss.secretKey));

            ZZ_pE ithPowerMultAlpha = ss.alpha * ithPower;
            powersOfSMultAlpha.append(encode(ithPowerMultAlpha, ss.secretKey));
            
            acc = acc * ss.s;
        }
    }


    long sizeOfImid = prog.circuit.numberOfMidWires;
    Vec<JLEncoding> rvVofS, rwWofS, ryYofS, alpharvVofS, alpharwWofS, alpharyYofS, betaSums;
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
        rvVofS[k] = encode(rvVkofS, ss.secretKey);
        rwWofS[k] = encode(rwWkofS, ss.secretKey);
        ryYofS[k] = encode(ryYkofS, ss.secretKey);
            
        //{E(alpha_v * r_v * v_k(S))}_k\in I_mid
        //{E(alpha_w * r_w * w_k(S))}_k\in I_mid
        //{E(alpha_y * r_y * y_k(S))}_k\in I_mid
        ZZ_pE alpharvVkofS, alpharwWkofS, alpharyYkofS;
        mul(alpharvVkofS, rvVkofS, ss.alpha_v);
        mul(alpharwWkofS, rwWkofS, ss.alpha_w);
        mul(alpharyYkofS, ryYkofS, ss.alpha_y);
        alpharvVofS[k] = encode(alpharvVkofS, ss.secretKey);
        alpharwWofS[k] = encode(alpharwWkofS, ss.secretKey);
        alpharyYofS[k] = encode(alpharyYkofS, ss.secretKey);

        //{E(beta ((r_v * v_k(S)) + (r_w * w_k(S)) + (r_y * y_k(S)))}_k\in I_mid
        ZZ_pE kthBetaSum;
        kthBetaSum = ss.beta * (rvVkofS + rwWkofS + ryYkofS);
        betaSums[k] = encode(kthBetaSum, ss.secretKey);
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
    crs.publicKey = ss.secretKey;//todo separate
    return crs;
}

Vec<ZZ> vecadd(const Vec<ZZ>& a, const Vec<ZZ>& b) {
    Vec<ZZ> res;
    if (a.length() < b.length()) {
        res = Vec<ZZ>(b);
        for (int i = 0; i < a.length(); i++) {
            res[i] += a[i];
        }
    } else {
        res = Vec<ZZ>(a);
        for (int i = 0; i < b.length(); i++) {
            res[i] += b[i];
        }
    }
    return res;
}

Proof prove(QRP prog, CRS crs, Vec<ZZ_p> allWireValues) {


    // Compute p = V*W-Y
    // P = W * W * Y = (Sum c_k * v_k(x)) * (Sum c_k * w_k(x)) - (Sum c_k * y_k(x))
    ZZ_pEX V, W, Y;
    for (int k = 0; k < prog.circuit.numberOfWires; k++) {
        V += allWireValues[k] * prog.V[k];
        W += allWireValues[k] * prog.W[k];
        Y += allWireValues[k] * prog.Y[k];
        ZZ_pEX Piter = V * W - Y;
        // if (IsZero(Piter)) cout << "0" << endl;
        // else  cout << "not 0" << endl;
    }
    // cout << V*W << "V*W" << endl;
    ZZ_pEX P = V*W-Y;
    // cout << P << "P" << endl;

    // Compute h = p / t
    ZZ_pEX H = P / prog.t;


    // E(r_v * Vmid(S))
    // E(r_v * Vmid(S) * alpha_v)
    // E(r_w * Wmid(S))
    // E(r_w * Wmid(S) * alpha_w)
    // E(r_y * Ymid(S))
    // E(r_y * Ymid(S) * alpha_y)
    // E(beta( (r_v * Vmid(S)) + (r_w * Wmid(S)) +(r_y * Ymid(S)) ))
    JLEncoding rvVmidOfS;
    JLEncoding rwWmidOfS;
    JLEncoding ryYmidOfS;
    JLEncoding alphaVrvVmidOfS;
    JLEncoding alphaWrwWmidOfS;
    JLEncoding alphaYryYmidOfS;
    JLEncoding betaSum;
    for (int k = prog.midOffset; k < prog.outOffset; k++) {
        //Scalar multiplications
        JLEncoding tmp;
        tmp = jle_scalar_mult(crs.rvVofS[k - prog.midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(rvVmidOfS, tmp, crs.publicKey);
        tmp = jle_scalar_mult(crs.rwWofS[k - prog.midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(rwWmidOfS, tmp, crs.publicKey);
        tmp = jle_scalar_mult(crs.ryYofS[k - prog.midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(ryYmidOfS, tmp, crs.publicKey);
        // rvVmidOfS += rep(allWireValues[k]) * crs.rvVofS[k - prog.midOffset];
        // rwWmidOfS += rep(allWireValues[k]) * crs.rwWofS[k - prog.midOffset];
        // ryYmidOfS += rep(allWireValues[k]) * crs.ryYofS[k - prog.midOffset];
        tmp = jle_scalar_mult(crs.alpharvVofS[k - prog.midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(alphaVrvVmidOfS, tmp, crs.publicKey);
        tmp = jle_scalar_mult(crs.alpharwWofS[k - prog.midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(alphaWrwWmidOfS, tmp, crs.publicKey);
        tmp = jle_scalar_mult(crs.alpharyYofS[k - prog.midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(alphaYryYmidOfS, tmp, crs.publicKey);
        // alphaVrvVmidOfS += rep(allWireValues[k]) * crs.alpharvVofS[k - prog.midOffset];
        // alphaWrwWmidOfS += rep(allWireValues[k]) * crs.alpharwWofS[k - prog.midOffset];
        // alphaYryYmidOfS += rep(allWireValues[k]) * crs.alpharyYofS[k - prog.midOffset];
        tmp = jle_scalar_mult(crs.betaSums[k - prog.midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(betaSum, tmp, crs.publicKey);
        // betaSum += rep(allWireValues[k]) * crs.betaSums[k - prog.midOffset];

    }

    // E(h(s))
    // E(alpha * h(s))

    // cout << H << "H" << endl;
    JLEncoding vec_hOfS, vec_alphaHofS;
    for (int i = 0 ; i <= deg(H); i++) {
        //todo How do we know that the degree of H is at most the number of multiplication gates?
        Vec<ZZ> ithCoeffOfH = to_vec_ZZ(rep(coeff(H, i)).rep);
        JLEncoding ihs = PlainMulEncryption(crs.powersOfS[i], ithCoeffOfH, crs.publicKey);
        JLEncoding iahs = PlainMulEncryption(crs.powersOfSMultAlpha[i], ithCoeffOfH, crs.publicKey);
        // cout << ihs << "ihs"<< endl;
        jle_add_assign(vec_hOfS, ihs, crs.publicKey);
        jle_add_assign(vec_alphaHofS, iahs, crs.publicKey);
        //todo should we reduce the degree of the encodings - and how?
    }

    return Proof{
        .rvVmidOfS = rvVmidOfS,
        .rwWmidOfS = rwWmidOfS,
        .ryYmidOfS  =ryYmidOfS,
        .alphaVrvVmidOfS = alphaVrvVmidOfS,
        .alphaWrwWmidOfS = alphaWrwWmidOfS,
        .alphaYryYmidOfS = alphaYryYmidOfS,
        .betaSum = betaSum,
        .hOfS = vec_hOfS,
        .alphaHOfS = vec_alphaHofS,
    };
}

bool verify(QRP qrp, secretState secret, CRS crs, Proof pi, Vec<ZZ_p> input, Vec<ZZ_p> output) {
    ZZ_pE rvVmidOfS = decode(pi.rvVmidOfS, secret.secretKey);
    ZZ_pE rwWmidOfS = decode(pi.rwWmidOfS, secret.secretKey);
    ZZ_pE ryYmidOfS = decode(pi.ryYmidOfS, secret.secretKey);
    ZZ_pE alphavrvVmidOfS = decode(pi.alphaVrvVmidOfS, secret.secretKey);
    ZZ_pE alphawrwWmidOfS = decode(pi.alphaWrwWmidOfS, secret.secretKey);
    ZZ_pE alphayryYmidOfS = decode(pi.alphaYryYmidOfS, secret.secretKey);
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
    ZZ_pE betaSum = decode(pi.betaSum, secret.secretKey);
    ZZ_pE computedBetaSum = secret.beta * (rvVmidOfS + rwWmidOfS + ryYmidOfS);

    if (computedBetaSum != betaSum) {
        cout << "beta sum \n";
        return false;
    }
    
    ZZ_pE hOfS = decode(pi.hOfS, secret.secretKey);
    ZZ_pE alphaHOfS = decode(pi.alphaHOfS, secret.secretKey);
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

    //r_y = r_v * r_w
    ZZ_pE computedV = secret.r_v * v_io + rvVmidOfS;
    ZZ_pE computedW = secret.r_w * w_io + rwWmidOfS;
    ZZ_pE computedY = secret.r_y * y_io + ryYmidOfS;

    ZZ_pE computedP = computedV * computedW - computedY;
    ZZ_pE hMultT = hOfS * eval(qrp.t, secret.s);
    // cout << "hOfS" << hOfS << endl;
    // cout << "hmult" << hMultT << endl;
    // cout << "ry" << secret.r_y << endl;
    if (computedP != secret.r_y * hMultT) {
        cout << "compute P: " << computedP << "\n";
        cout << "r_y * h*t: " << secret.r_y * hMultT << "\n";
        return false;
    }

    return  true;
}
