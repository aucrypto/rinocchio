#include <NTL/ZZ.h>

using namespace NTL;

void keygen(ZZ& p_prime, ZZ& p, ZZ& n, ZZ& y, ZZ& D, const int l, const int k);

void precompute_pow2(ZZ& pow2k, const int k);

void encrypt(ZZ& c, const ZZ& m, const ZZ& n, const ZZ& y, const int k, const ZZ& pow2k);

void decrypt(ZZ& m, const ZZ& c, const ZZ& p, const ZZ& p_prime, const ZZ& D_, const int k, const ZZ& pow2k1);

void scalar_mult_encrypted(ZZ& c_result, const ZZ& c, const ZZ& scalar, ZZ& n);

void add_encrypted(ZZ& c_result, ZZ& c1, const ZZ& c2, const ZZ& n);