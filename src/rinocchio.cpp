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

Proof prove(const QRP& prog, const CRS& crs, const Vec<ZZ_p>& allWireValues) {
    cout << "Start prove...\n";
    clock_t t;
    // Compute p = V*W-Y
    // P = W * W * Y = (Sum c_k * v_k(x)) * (Sum c_k * w_k(x)) - (Sum c_k * y_k(x))
    t = clock();
    ZZ_pEX V, W, Y;
    for (int k = 0; k < prog.circuit.numberOfWires; k++) {
        V += allWireValues[k] * prog.V[k];
        W += allWireValues[k] * prog.W[k];
        Y += allWireValues[k] * prog.Y[k];
    }
    t = clock() - t;
    cout << "V, W and Y computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
    // cout << V*W << "V*W" << endl;

    t = clock();
    ZZ_pEX P = V*W-Y;
    // cout << P << "P" << endl;

    // Compute h = p / t
    ZZ_pEX H = P / prog.t;
    t = clock() - t;
    cout << "H computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";


    // E(r_v * Vmid(S))
    // E(r_v * Vmid(S) * alpha_v)
    // E(r_w * Wmid(S))
    // E(r_w * Wmid(S) * alpha_w)
    // E(r_y * Ymid(S))
    // E(r_y * Ymid(S) * alpha_y)
    // E(beta( (r_v * Vmid(S)) + (r_w * Wmid(S)) +(r_y * Ymid(S)) ))
    t = clock();
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
    t = clock() - t;
    cout << "mid polynomials computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    // E(h(s))
    // E(alpha * h(s))

    // cout << H << "H" << endl;
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
        jle_add_assign(vec_hOfS, ihs, crs.publicKey);
        jle_add_assign(vec_alphaHofS, iahs, crs.publicKey);
    }
    t = clock() - t;
    cout << "H(s) computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";//slow

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

bool verify(const SecretState& secret, const CRS& crs, const Proof& pi, const Vec<ZZ_p>& input, const Vec<ZZ_p>& output) {
    cout << "Verify start\n";
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
    cout << "Compute io polys\n";
    t = clock();
    // compute P_in
    ZZ_pE v_io, w_io, y_io;
    for (int k = 0; k < input.length(); k++) {
        v_io += input[k] * secret.inVofS[k];
        w_io += input[k] * secret.inWofS[k];
        y_io += input[k] * secret.inYofS[k];
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
    ZZ_pE computedV = secret.r_v * v_io + rvVmidOfS;
    ZZ_pE computedW = secret.r_w * w_io + rwWmidOfS;
    ZZ_pE computedY = secret.r_y * y_io + ryYmidOfS;

    ZZ_pE computedP = computedV * computedW - computedY;
    ZZ_pE hMultT = hOfS * secret.tOfS;// todo eliminate eval
    t = clock() - t;
    cout << "done checking h, took " << ((double) t) / CLOCKS_PER_SEC << "\n";
    if (computedP != secret.r_y * hMultT) {
        cout << "compute P: " << computedP << "\n";
        cout << "r_y * h*t: " << secret.r_y * hMultT << "\n";
        return false;
    }
    cout << "Verified!\n";
    return  true;
}
