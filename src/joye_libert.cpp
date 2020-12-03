#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <vector>
#include <joye_libert.h>

using namespace std;
using namespace NTL;

void keygen(ZZ& p_prime, ZZ& p, ZZ& n, ZZ& g, ZZ& D, const int l, const int k) {
    ZZ p_;
    bool p_is_prime = false;
    while (!p_is_prime) {
        p_prime = GenPrime_ZZ(l-k, 10);
        p = p_prime << k;
        p = p + 1;
        long primality_test = ProbPrime(p);
        // cout << "candidate: " << p << "\n";
        // cout << "primality test: " << primality_test << "\n";

        // TODO: make sure that p is prime with high enough probability.
        // maybe increase number of Miller-Rabin witness tests.
        p_is_prime = (primality_test == 1); 
    }

    ZZ q = GenGermainPrime_ZZ(l-1);
    q = q << 1;
    q = q + 1;

    n = p * q;

    bool y_is_gen = false;
    while (!y_is_gen) {
        g = RandomBnd(n);
        bool is_gen_in_Zp = Jacobi(g, p) == -1;
        bool is_gen_in_Zq = Jacobi(g, q) == -1;

        y_is_gen = is_gen_in_Zp && is_gen_in_Zq;
    }

    ZZ y_mod_p = g % p;
    D = PowerMod(y_mod_p, p_prime, p);
    D = InvMod(D, p);
}

/*
 * Computes and assigns 2^k to pow2k
 */
void precompute_pow2(ZZ& pow2k, const int k) {
    pow2k = ZZ(1) << k;
}

/*
 * Encrypts input message m and assigns the resulting cipher to c
 * Parameters:
 *   - n, g, k: public parameters (n, g, k)
 *   - pow2k: precomputed 2^k
 *   - c: pointer to which output is assigned
 * Precondition:
 *   - n, g, k, pow2k are valid key and power of 2.
 */
void encrypt(ZZ& c, const ZZ& m, const ZZ& n, const ZZ& g, const int k, const ZZ& pow2k) {
	ZZ x;
    x = RandomBnd(n);
    x = PowerMod(x, pow2k, n);
    c = PowerMod(g, m, n);
    c = (c * x) % n;
}

/*
 * Decrypts cipher c and assigns the resulting message to m
 * Parameters:
 *   - p_prime: private key
 *   - p: precomputed p*2^k + 1
 *   - D_: precomputed g^(-p_prime)
 *   - k: bit-length of m
 *   - pow2k1: precomputed 2^(k-1)
 * Precondition:
 *   - p_prime, p, D_, k are all computed correctly
 */
void decrypt(ZZ& m, const ZZ& c, const ZZ& p, const ZZ& p_prime, const ZZ& D_, const int k, const ZZ& pow2k1) {
    ZZ C = c % p;
    C = PowerMod(C, p_prime, p);
    m = ZZ(0);
    
    ZZ B, D, E, temp;
    int i;
    B = ZZ(1);
    D = D_;
    E = pow2k1;
    for (i=0;i<k-1;i++) {
        temp = PowerMod(C, E, p);
        if (IsOne(temp) != 1) {
            m = m + B;
            temp = C * D;
            C = temp % p;
        }
        B = B + B;
        D = PowerMod(D, 2, p);
        E = E >> 1;
    }

    if (IsOne(C) != 1) {
        m = m + B;
    }
}

void add_encrypted(ZZ& result, const ZZ& c1, const ZZ& c2, const ZZ& n) {
    result = c1 * c2;
    result = result % n;
}

void scalar_mult_encrypted(ZZ& c_result, const ZZ& c, const ZZ& scalar, const ZZ& n) {
    c_result = PowerMod(c, scalar, n);
}

JLEncodingKey gen_jl_encoding_key(long l, long k) {
    ZZ p_prime, p, n, g, D;
    keygen(p_prime, p, n, g, D, l, k);

    ZZ pow2k, pow2k1;
    precompute_pow2(pow2k, k);
    precompute_pow2(pow2k1, k-1);

    return JLEncodingKey{
        .p_prime = p_prime,
        .p = p,
        .n = n,
        .g = g,
        .D = D,
        .pow2k = pow2k, 
        .pow2k1 = pow2k1,
        .l = l, 
        .k = k, 
    };
}

JLEncoding encode(const ZZ_pE& m, const JLEncodingKey& key) {
    Vec<ZZ> res;
    const ZZ_pX& polynomialRep = rep(m);
    
    for (int i = 0; i < ZZ_pE::degree(); i++) {
        ZZ c;
        const ZZ_p& ci = coeff(polynomialRep, i);
        encrypt(c, rep(ci), key.n, key.g, key.k, key.pow2k);
        res.append(c);
    }
    return JLEncoding{.coeffs = res};
}

ZZ_pE decode(const JLEncoding& c, const JLEncodingKey& key) {
    ZZ_pX result;
    
    // std::cout << "k....: " << key.k << "\n";
    // std::cout << "po2k1: " << key.pow2k1 << "\n";
    ZZ ci;
    for (int i = 0; i < c.coeffs.length(); i++) {
        // std::cout << "coeff: " << c.coeffs[i] << "\n";
        // std::cout << "p....: " << key.p << "\n";
        // std::cout << "p_pri: " << key.p_prime << "\n";
        // std::cout << "D....: " << key.D << "\n";
        decrypt(ci, c.coeffs[i], key.p, key.p_prime, key.D, key.k, key.pow2k1);
        // std::cout << "decry: " << ci << "\n";
        ZZ_p ci_ = conv<ZZ_p>(ci);
        SetCoeff(result, i, ci_);
    }

    return conv<ZZ_pE>(result);
}

void jle_add_assign(JLEncoding& a, const JLEncoding& b, const JLEncodingKey& key) {
    for (int i = 0; i < a.coeffs.length() && i < b.coeffs.length(); i++) {
        add_encrypted(a.coeffs[i], a.coeffs[i], b.coeffs[i], key.n);
    }
    if (a.coeffs.length() < b.coeffs.length()) {
        int startIndex = a.coeffs.length();
        a.coeffs.SetLength(b.coeffs.length());
        for (int i = startIndex; i < a.coeffs.length(); i++) {
            a.coeffs[i] = b.coeffs[i];
        }
    }
}

void jle_scalar_mult_assign(JLEncoding& a, const ZZ& scalar, const JLEncodingKey& key) {
    for (int i = 0; i < a.coeffs.length() ; i++) {
        scalar_mult_encrypted(a.coeffs[i], a.coeffs[i], scalar, key.n);
    }
}

JLEncoding jle_mult(const JLEncoding& a, const Vec<ZZ>& b, const JLEncodingKey& key) {
    Vec<ZZ> res;
    int d = a.coeffs.length() + b.length() - 1;
    res.SetLength(d); //todo -2?
    for (int i = 0; i < d; i++) {
        res[i] = ZZ(1);
    }
    for (int i = 0; i < a.coeffs.length(); i++) {
        for (int j = 0; j < b.length(); j++) {
            ZZ tmp;
            scalar_mult_encrypted(tmp, a.coeffs[i], b[j], key.n);
            add_encrypted(res[i+j], res[i+j], tmp, key.n);
        }
    }
    return JLEncoding{.coeffs = res}; 
}
//todo compare multiplications
JLEncoding PlainMulEncryption(const JLEncoding& a, const Vec<ZZ>& b, JLEncodingKey& key) {
    JLEncoding res;
    long da = a.coeffs.length() - 1;
    long db = b.length() - 1;
    long d = da+db;
    res.coeffs.SetLength(d + 1); //todo lengths - 1



    const ZZ *ap, *bp;

    ap = a.coeffs.elts();
    bp = b.elts();

    ZZ *resp = res.coeffs.elts();

    long i, j, jmin, jmax;
    ZZ t, acc;

    for (i = 0; i <= d; i++) {
        jmin = max(0, i-db);
        jmax = min(da, i);
        set(acc);
        for (j = jmin; j <= jmax; j++) {
            scalar_mult_encrypted(t, ap[j], bp[i-j], key.n);
            add_encrypted(acc, acc, t, key.n);
        }
        resp[i] = acc;
    }
    return res;
}

ostream& operator<<(ostream& s, const JLEncoding& jle) {
    return s << jle.coeffs;
}
