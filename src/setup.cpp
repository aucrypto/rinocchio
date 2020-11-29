#include <setup.h>
#include <gr.h>
#include <joye_libert.h>
#include <NTL/ZZ_pE.h>

using namespace NTL;

secretState setup(long l, long k) {
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
    ss.secretKey = gen_jl_encoding_key(l, k);
    
    return ss;
}