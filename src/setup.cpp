#include <setup.h>
#include <gr.h>
#include <joye_libert.h>
#include <NTL/ZZ_pE.h>

using namespace NTL;

SecretState setup(const QRP& qrp, long l, long k) {
    SecretState ss;

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

    //todo we could just as well scale by r_v, r_w and r_y in the precomputations, would save 300 ZZ_pE mults in 10 dimension case when verifying
    ss.inVofS.SetLength(qrp.circuit.numberOfInputWires);
    ss.inWofS.SetLength(qrp.circuit.numberOfInputWires);
    ss.inYofS.SetLength(qrp.circuit.numberOfInputWires);
    for (int k = 0; k < qrp.circuit.numberOfInputWires; k++) {
        ss.inVofS[k] = eval(qrp.V[k], s);
        ss.inWofS[k] = eval(qrp.W[k], s);
        ss.inYofS[k] = eval(qrp.Y[k], s);
    }
    ss.outVofS.SetLength(qrp.circuit.numberOfOutputWires);
    ss.outWofS.SetLength(qrp.circuit.numberOfOutputWires);
    ss.outYofS.SetLength(qrp.circuit.numberOfOutputWires);
    for (int k = 0; k < qrp.circuit.numberOfOutputWires; k++) {
        ss.outVofS[k] = eval(qrp.V[k + qrp.outOffset], s);
        ss.outWofS[k] = eval(qrp.W[k + qrp.outOffset], s);
        ss.outYofS[k] = eval(qrp.Y[k + qrp.outOffset], s);
    }

    ss.s = s;
    ss.r_v = r_v;
    ss.r_w = r_w;
    ss.r_y = r_y;
    ss.alpha = alpha;
    ss.alpha_v = alpha_v;
    ss.alpha_w = alpha_w;
    ss.alpha_y = alpha_y;
    ss.beta = beta;
    ss.tOfS = eval(qrp.t, s);
    ss.secretKey = gen_jl_encoding_key(l, k);

    return ss;
}


CRS getCRS(const QRP& prog, const SecretState& ss) {
    
    //{E(S^i)}_i=0^d
    //{E(alpha * S^i)}_i=0^d
    Vec<JLEncoding> powersOfS;
    Vec<JLEncoding> powersOfSMultAlpha;
    {
        ZZ_pE acc;
        set(acc);
        for (int i = 0; i <= prog.circuit.numberOfMultiplicationGates; i++) {
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

    CRS crs;//Todo is it significantly faster not to copy here?
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


ostream& operator<<(ostream& s, const SecretState& secretState)  {
    return s;
}

istream& operator>>(istream& s, SecretState& secretState) {
    return s;
}


ostream& operator<<(ostream& s, const CRS& crs) {
    return s;
}

istream& operator>>(istream& s, CRS& crs) {
    return s;
}
