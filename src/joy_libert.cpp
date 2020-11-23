#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

void keygen(ZZ& p, ZZ& n, ZZ& y, ZZ& D, const int l, const int k) {
    ZZ p_;
    bool p__is_prime = false;
    while (!p__is_prime) {
        p = GenPrime_ZZ(l-k, 10);
        p_ = p << k;
        p_ = p_ + 1;
        long primality_test = ProbPrime(p_);
        // cout << "candidate: " << p_ << "\n";
        // cout << "primality test: " << primality_test << "\n";
        p__is_prime = (primality_test == 1); // TODO: make sure that p_ is prime with high enough probability.
                                            // maybe increase number of miller-rabin witness tests.
    }
    cout << "found p': " << p << "\n";
    cout << "found p: " << p_ << "\n";

    ZZ q_ = GenGermainPrime_ZZ(l-1);
    q_ = q_ << 1;
    q_ = q_ + 1;
    cout << "found q: " << q_ << "\n";

    n = p_ * q_;
    cout << "found n: " << n << "\n";

    bool y_is_gen = false;
    while (!y_is_gen) {
        y = RandomBnd(n);
        bool is_gen_in_Zp = Jacobi(y, p_) == -1;
        bool is_gen_in_Zq = Jacobi(y, q_) == -1;

        y_is_gen = is_gen_in_Zp && is_gen_in_Zq;
    }
    cout << "found g: " << y << "\n";

    ZZ y_mod_p = y % p;
    D = PowerMod(y_mod_p, p, p_);
    D = InvMod(D, p_);
    cout << "found D: " << D << "\n";
}

// int precompute(mpz_t _2k1, mpz_t _2k, mpz_t pm12k, 
// 	         const mpz_t p, const int k);
