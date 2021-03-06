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

#include <time.h>

using namespace std;
using namespace NTL;

ZZ_pEX proverComputeH(const QRP& prog, const Vec<ZZ_p>& allWireValues) {
    ZZ_pEX H;
    clock_t t;
    // Compute p = V*W-Y
    // P = W * W * Y = (Sum c_k * v_k(x)) * (Sum c_k * w_k(x)) - (Sum c_k * y_k(x))
    t = clock();
    ZZ_pEX V, W, Y;
    for (int k = 0; k < prog.circuit.numberOfInputWires; k++) {
        V += allWireValues[k] * prog.V[k];
        W += allWireValues[k] * prog.W[k];
    }

    for (int k = prog.circuit.numberOfInputWires; k < prog.circuit.numberOfWires; k++) {
        V += allWireValues[k] * prog.V[k];
        W += allWireValues[k] * prog.W[k];
        Y += allWireValues[k] * prog.Y[k];//todo handle input = 0
    }
    t = clock() - t;
    cout << "V, W and Y computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    
    t = clock();
    ZZ_pEX P = V*W-Y;
    // cout << P << "P" << endl;

    // Compute h = p / t
    H = P / prog.t;
    t = clock() - t;
    cout << "H computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
    return H;
}

Proof prove(const ZZ_pEX& H, const CRS& crs, const Vec<ZZ_p>& allWireValues, long midOffset, long outOffset) {
    Proof proof;
    clock_t t;

    // E(r_v * Vmid(S))
    // E(r_v * Vmid(S) * alpha_v)
    // E(r_w * Wmid(S))
    // E(r_w * Wmid(S) * alpha_w)
    // E(r_y * Ymid(S))
    // E(r_y * Ymid(S) * alpha_y)
    // E(beta( (r_v * Vmid(S)) + (r_w * Wmid(S)) +(r_y * Ymid(S)) ))
    t = clock();
    for (int k = midOffset; k < outOffset; k++) {//todo entire loop could use length of rvVofS instead of referencing qrp.
        //Scalar multiplications
        JLEncoding tmp;
        tmp = jle_scalar_mult(crs.rvVofS[k - midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(proof.rvVmidOfS, tmp, crs.publicKey);
        tmp = jle_scalar_mult(crs.rwWofS[k - midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(proof.rwWmidOfS, tmp, crs.publicKey);
        tmp = jle_scalar_mult(crs.ryYofS[k - midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(proof.ryYmidOfS, tmp, crs.publicKey);
        // rvVmidOfS += rep(allWireValues[k]) * crs.rvVofS[k - prog.midOffset];
        // rwWmidOfS += rep(allWireValues[k]) * crs.rwWofS[k - prog.midOffset];
        // ryYmidOfS += rep(allWireValues[k]) * crs.ryYofS[k - prog.midOffset];
        tmp = jle_scalar_mult(crs.alpharvVofS[k - midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(proof.alphaVrvVmidOfS, tmp, crs.publicKey);
        tmp = jle_scalar_mult(crs.alpharwWofS[k - midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(proof.alphaWrwWmidOfS, tmp, crs.publicKey);
        tmp = jle_scalar_mult(crs.alpharyYofS[k - midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(proof.alphaYryYmidOfS, tmp, crs.publicKey);
        // alphaVrvVmidOfS += rep(allWireValues[k]) * crs.alpharvVofS[k - prog.midOffset];
        // alphaWrwWmidOfS += rep(allWireValues[k]) * crs.alpharwWofS[k - prog.midOffset];
        // alphaYryYmidOfS += rep(allWireValues[k]) * crs.alpharyYofS[k - prog.midOffset];
        tmp = jle_scalar_mult(crs.betaSums[k - midOffset], rep(allWireValues[k]), crs.publicKey);
        jle_add_assign(proof.betaSum, tmp, crs.publicKey);
        // betaSum += rep(allWireValues[k]) * crs.betaSums[k - prog.midOffset];

    }
    t = clock() - t;
    cout << "Encryptions of mid wire polynomials computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    // E(h(s))
    // E(alpha * h(s))

    t = clock();
    JLEncoding vec_hOfS, vec_alphaHofS;
    for (int i = 0 ; i <= deg(H); i++) {
        //todo How do we know that the degree of H is at most the number of multiplication gates?
        Vec<ZZ> ithCoeffOfH = to_vec_ZZ(rep(coeff(H, i)).rep); //todo is this often zero?
        // cout << i << "th coeff of H:" << ithCoeffOfH << endl;
        JLEncoding ihs = PlainMulEncryption(crs.powersOfS[i], ithCoeffOfH, crs.publicKey.n);
        JLEncoding iahs = PlainMulEncryption(crs.powersOfSMultAlpha[i], ithCoeffOfH, crs.publicKey.n);
        //Todo computing the reduction would make the next step twice as fast. (And a shorter proof)
        // cout << ihs << "ihs"<< endl;
        jle_add_assign(proof.hOfS, ihs, crs.publicKey);
        jle_add_assign(proof.alphaHOfS, iahs, crs.publicKey);
    }
    t = clock() - t;
    cout << "H(s) computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";//slow

    return proof;
}

Proof prove(const QRP& qrp, const CRS& crs, const Vec<ZZ_p>& allWireValues) {
    // Proof proof;
    clock_t t;
    ZZ_pEX H = proverComputeH(qrp, allWireValues);
    t = clock() - t;
    cout << "H computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    return prove(H, crs, allWireValues, qrp.midOffset, qrp.outOffset);
}

bool verify(const SecretState& secret, const CRS& crs, const Proof& pi, const Vec<ZZ_p>& input, const Vec<ZZ_p>& output) {
    clock_t t = clock();
    ZZ_pE rvVmidOfS = decode(pi.rvVmidOfS, secret.secretKey);
    ZZ_pE rwWmidOfS = decode(pi.rwWmidOfS, secret.secretKey);
    ZZ_pE ryYmidOfS = decode(pi.ryYmidOfS, secret.secretKey);
    ZZ_pE alphavrvVmidOfS = decode(pi.alphaVrvVmidOfS, secret.secretKey);
    ZZ_pE alphawrwWmidOfS = decode(pi.alphaWrwWmidOfS, secret.secretKey);
    ZZ_pE alphayryYmidOfS = decode(pi.alphaYryYmidOfS, secret.secretKey);
    ZZ_pE betaSum = decode(pi.betaSum, secret.secretKey);
    ZZ_pE hOfS = decode(pi.hOfS, secret.secretKey);
    ZZ_pE alphaHOfS = decode(pi.alphaHOfS, secret.secretKey);
    t = clock() - t;
    cout << "done decoding, took " << ((double) t) / CLOCKS_PER_SEC << "\n";
    t = clock();
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
    ZZ_pE computedBetaSum = secret.beta * (rvVmidOfS + rwWmidOfS + ryYmidOfS);

    if (computedBetaSum != betaSum) {
        cout << "beta sum \n";
        return false;
    }
    
    if (secret.alpha * hOfS != alphaHOfS) {
        cout << "\n";
        return false;
    }

    t = clock() - t;
    cout << "done checking decodings, took " << ((double) t) / CLOCKS_PER_SEC << "\n";
    t = clock();
    // compute P_in
    ZZ_pE v_io, w_io, y_io;
    for (int k = 0; k < input.length(); k++) {
        v_io += input[k] * secret.inVofS[k];
        w_io += input[k] * secret.inWofS[k];
        //The y polynomial for inputs is always zero
    }
    // add P_out
    for (int k = 0; k < output.length(); k++) {
        v_io += output[k] * secret.outVofS[k];
        w_io += output[k] * secret.outWofS[k];
        y_io += output[k] * secret.outYofS[k];
    }
    t = clock() - t;
    cout << "done computing io, took " << ((double) t) / CLOCKS_PER_SEC << "\n";
    t = clock();
    //r_y = r_v * r_w
    ZZ_pE computedV = v_io + rvVmidOfS;
    ZZ_pE computedW = w_io + rwWmidOfS;
    ZZ_pE computedY = y_io + ryYmidOfS;

    ZZ_pE computedP = computedV * computedW - computedY;
    ZZ_pE hMultT = hOfS * secret.tOfS;
    t = clock() - t;
    cout << "done checking h, took " << ((double) t) / CLOCKS_PER_SEC << "\n";
    if (computedP != hMultT) {
        cout << "compute P: " << computedP << "\n";
        cout << "r_y * h*t: " << secret.r_y * hMultT << "\n";
        return false;
    }
    return  true;
}
