#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

void keygen(ZZ& p_prime, ZZ& p, ZZ& n, ZZ& y, ZZ& D, const int l, const int k) {
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
        y = RandomBnd(n);
        bool is_gen_in_Zp = Jacobi(y, p) == -1;
        bool is_gen_in_Zq = Jacobi(y, q) == -1;

        y_is_gen = is_gen_in_Zp && is_gen_in_Zq;
    }

    ZZ y_mod_p = y % p;
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
 *   - n, y, k: public parameters (n, g, k)
 *   - pow2k: precomputed 2^k
 *   - c: pointer to which output is assigned
 * Precondition:
 *   - n, y, k, pow2k are valid key and power of 2.
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

void add_encrypted(ZZ& c_result, ZZ& c1, ZZ& c2, const ZZ& n) {
    c_result = c1 * c2;
    c_result = c_result % n;
}

void scalar_mult_encrypted(ZZ& c_result, ZZ& c, ZZ& scalar, ZZ& n) {
    c_result = PowerMod(c, scalar, n);
}