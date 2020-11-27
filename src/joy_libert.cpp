#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <vector>

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

struct JLEncodingKeyPart {
    ZZ p_prime, p, n, g, D;
};

struct JLEncodingKey {
    Vec<JLEncodingKeyPart> keys;
    ZZ pow2k, pow2k1;
    long l, k, d;
};

JLEncodingKey gen_jl_encoding_key(long l, long k, long d) {
    Vec<JLEncodingKeyPart> keys;
    for (int i = 0; i < d; i++) {
        ZZ p_prime, p, n, g, D;
        keygen(p_prime, p, n, g, D, l, k);
        keys.append(JLEncodingKeyPart{
            p_prime: p_prime,
            p: p,
            n: n,
            g: g,
            D: D,
        });
    };

    ZZ pow2k, pow2k1;
    precompute_pow2(pow2k, k);
    precompute_pow2(pow2k1, k-1);

    return JLEncodingKey{
        keys: keys,
        pow2k: pow2k, 
        pow2k1: pow2k1,
        l: l, 
        k: k, 
        d: d,
    };
}


struct JLEncoding {
    Vec<ZZ> coeffs;
};

JLEncoding encode(const ZZ_pE& m, const JLEncodingKey& key) {
    Vec<ZZ> res;
    const ZZ_pX& polynomialRep = rep(m);
    
    for (int i = 0; i < ZZ_pE::degree(); i++) {
        ZZ c;
        const ZZ_p& ci = polynomialRep[i];
        const JLEncodingKeyPart& keypart_i = key.keys[i];
        encrypt(c, rep(ci), keypart_i.n, keypart_i.g, key.k, key.pow2k);
        res.append(c);
    }
    return JLEncoding{coeffs: res};
}

ZZ_pE decode(const JLEncoding& c, const JLEncodingKey& key) {
    ZZ_pX result;
    
    std::cout << "k....: " << key.k << "\n";
    std::cout << "po2k1: " << key.pow2k1 << "\n";
    ZZ ci;
    for (int i = 0; i < key.d; i++) {
        const JLEncodingKeyPart& keypart_i = key.keys[i];
        std::cout << "coeff: " << c.coeffs[i] << "\n";
        std::cout << "p....: " << keypart_i.p << "\n";
        std::cout << "p_pri: " << keypart_i.p_prime << "\n";
        std::cout << "D....: " << keypart_i.D << "\n";
        decrypt(ci, c.coeffs[i], keypart_i.p, keypart_i.p_prime, keypart_i.D, key.k, key.pow2k1);
        std::cout << "decry: " << ci << "\n";
        ZZ_p ci_ = conv<ZZ_p>(ci);
        SetCoeff(result, i, ci_);
    }

    return conv<ZZ_pE>(result);
}

void jle_add_assign(JLEncoding& a, const JLEncoding& b, const JLEncodingKey& key) {
    for (int i = 0; i < key.d; i++) {
        add_encrypted(a.coeffs[i], a.coeffs[i], b.coeffs[i], key.keys[i].n);
    }
}

void jle_mult_assign(JLEncoding& a, const ZZ& scalar, const JLEncodingKey& key) {
    for (int i = 0; i < key.d; i++) {
        scalar_mult_encrypted(a.coeffs[i], a.coeffs[i], scalar, key.keys[i].n);
    }
}