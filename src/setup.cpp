#include <setup.h>
#include <gr.h>
#include <joye_libert.h>
#include <NTL/ZZ_pE.h>

#include <time.h>

using namespace NTL;

SecretState setup(const QRP& qrp, long l, long k) {
    SecretState ss;
    clock_t t = clock();
    //Find non-zero s in exceptional set (i.e. in A^*):
    ss.s = randomNonZeroInExceptionalSet();
    //Find random r_v, r_w in R^*, i.e. d random coefficients where at least one is odd
    ss.r_v = randomInvertible();
    ss.r_w = randomInvertible();
    ss.r_y = ss.r_v * ss.r_w;
    
    // ZZ_pE alpha, alpha_v, alpha_w, alpha_y;
    ss.alpha = randomInvertible();
    ss.alpha_v = randomInvertible();
    ss.alpha_w = randomInvertible();
    ss.alpha_y = randomInvertible();

    do {
        ss.beta = random_ZZ_pE();
    } while (IsZero(ss.beta));
    t = clock() - t;
    std::cout << "Secret elements found: " << ((double) t) / CLOCKS_PER_SEC << " seconds.\n";
    t = clock();

    ss.inVofS.SetLength(qrp.circuit.numberOfInputWires);
    ss.inWofS.SetLength(qrp.circuit.numberOfInputWires);
    for (int k = 0; k < qrp.circuit.numberOfInputWires; k++) {
        ss.inVofS[k] = eval(qrp.V[k], ss.s) * ss.r_v;
        ss.inWofS[k] = eval(qrp.W[k], ss.s) * ss.r_w;
    }
    ss.outVofS.SetLength(qrp.circuit.numberOfOutputWires);
    ss.outWofS.SetLength(qrp.circuit.numberOfOutputWires);
    ss.outYofS.SetLength(qrp.circuit.numberOfOutputWires);
    for (int k = 0; k < qrp.circuit.numberOfOutputWires; k++) {
        ss.outVofS[k] = eval(qrp.V[k + qrp.outOffset], ss.s) * ss.r_v;
        ss.outWofS[k] = eval(qrp.W[k + qrp.outOffset], ss.s) * ss.r_w;
        ss.outYofS[k] = eval(qrp.Y[k + qrp.outOffset], ss.s) * ss.r_y;
    }
    t = clock() - t;
    std::cout << "Wire polynomials evaluated in s: " << ((double) t) / CLOCKS_PER_SEC << " seconds.\n";


    ss.tOfS = eval(qrp.t, ss.s) * ss.r_y;
    ss.secretKey = gen_jl_encoding_key(l, k);

    return ss;
}


CRS getCRS(const QRP& prog, const SecretState& ss) {
    CRS crs;
    crs.publicKey = ss.secretKey;//todo separate pk
    clock_t t = clock();
    //{E(S^i)}_i=0^d
    //{E(alpha * S^i)}_i=0^d
    {
        ZZ_pE acc;
        set(acc);
        for (int i = 0; i <= prog.circuit.numberOfMultiplicationGates; i++) {
            ZZ_pE ithPower = ZZ_pE(acc);
            crs.powersOfS.append(encode(ithPower, ss.secretKey));

            ZZ_pE ithPowerMultAlpha = ss.alpha * ithPower;
            crs.powersOfSMultAlpha.append(encode(ithPowerMultAlpha, ss.secretKey));
            
            acc = acc * ss.s;
        }
    }

    t = clock() - t;
    std::cout << "Powers of s encrypted: " << ((double) t) / CLOCKS_PER_SEC << " seconds.\n";
    t = clock();

    long sizeOfImid = prog.circuit.numberOfMidWires;
    crs.rvVofS.SetLength(sizeOfImid);
    crs.rwWofS.SetLength(sizeOfImid);
    crs.ryYofS.SetLength(sizeOfImid);
    crs.alpharvVofS.SetLength(sizeOfImid);
    crs.alpharwWofS.SetLength(sizeOfImid);
    crs.alpharyYofS.SetLength(sizeOfImid);
    crs.betaSums.SetLength(sizeOfImid);
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
        crs.rvVofS[k] = encode(rvVkofS, ss.secretKey);
        crs.rwWofS[k] = encode(rwWkofS, ss.secretKey);
        crs.ryYofS[k] = encode(ryYkofS, ss.secretKey);
            
        //{E(alpha_v * r_v * v_k(S))}_k\in I_mid
        //{E(alpha_w * r_w * w_k(S))}_k\in I_mid
        //{E(alpha_y * r_y * y_k(S))}_k\in I_mid
        ZZ_pE alpharvVkofS, alpharwWkofS, alpharyYkofS;
        mul(alpharvVkofS, rvVkofS, ss.alpha_v);
        mul(alpharwWkofS, rwWkofS, ss.alpha_w);
        mul(alpharyYkofS, ryYkofS, ss.alpha_y);
        crs.alpharvVofS[k] = encode(alpharvVkofS, ss.secretKey);
        crs.alpharwWofS[k] = encode(alpharwWkofS, ss.secretKey);
        crs.alpharyYofS[k] = encode(alpharyYkofS, ss.secretKey);

        //{E(beta ((r_v * v_k(S)) + (r_w * w_k(S)) + (r_y * y_k(S)))}_k\in I_mid
        ZZ_pE kthBetaSum;
        kthBetaSum = ss.beta * (rvVkofS + rwWkofS + ryYkofS);
        crs.betaSums[k] = encode(kthBetaSum, ss.secretKey);
    }
    t = clock() - t;
    std::cout << "Masked wire polynomials in s encrypted: " << ((double) t) / CLOCKS_PER_SEC << " seconds.\n";

    return crs;
}


ostream& operator<<(ostream& s, const SecretState& secretState)  {
    s << secretState.s;
    s << "\n";
    s << secretState.r_v;
    s << "\n";
    s << secretState.r_w;
    s << "\n";
    s << secretState.r_y;
    s << "\n";
    s << secretState.alpha;
    s << "\n";
    s << secretState.alpha_v;
    s << "\n";
    s << secretState.alpha_w;
    s << "\n";
    s << secretState.alpha_y;
    s << "\n";
    s << secretState.beta;
    s << "\n";
    s << secretState.tOfS;
    s << "\n";
    s << secretState.inVofS;
    s << "\n";
    s << secretState.inWofS;
    s << "\n";
    s << secretState.outVofS;
    s << "\n";
    s << secretState.outWofS;
    s << "\n";
    s << secretState.outYofS;
    s << "\n";
    s << secretState.secretKey;
    s << "\n";
    return s;
}

istream& operator>>(istream& s, SecretState& secretState) {
    s >> secretState.s;
    s >> secretState.r_v;
    s >> secretState.r_w;
    s >> secretState.r_y;
    s >> secretState.alpha;
    s >> secretState.alpha_v;
    s >> secretState.alpha_w;
    s >> secretState.alpha_y;
    s >> secretState.beta;
    s >> secretState.tOfS;
    s >> secretState.inVofS;
    s >> secretState.inWofS;
    s >> secretState.outVofS;
    s >> secretState.outWofS;
    s >> secretState.outYofS;
    s >> secretState.secretKey;
    return s;
}


ostream& operator<<(ostream& s, const CRS& crs) {
    s << crs.powersOfS;
    s << "\n";
    s << crs.powersOfSMultAlpha;
    s << "\n";
    s << crs.rvVofS;
    s << "\n";
    s << crs.rwWofS;
    s << "\n";
    s << crs.ryYofS;
    s << "\n";
    s << crs.alpharvVofS;
    s << "\n";
    s << crs.alpharwWofS;
    s << "\n";
    s << crs.alpharyYofS;
    s << "\n";
    s << crs.betaSums; 
    s << "\n";
    s << crs.publicKey;
    s << "\n";

    return s;
}

istream& operator>>(istream& s, CRS& crs) {
    s >> crs.powersOfS;
    s >> crs.powersOfSMultAlpha;
    s >> crs.rvVofS;
    s >> crs.rwWofS;
    s >> crs.ryYofS;
    s >> crs.alpharvVofS;
    s >> crs.alpharwWofS;
    s >> crs.alpharyYofS;
    s >> crs.betaSums; 
    s >> crs.publicKey;
    return s;
}
