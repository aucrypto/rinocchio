#include <NTL/ZZ.h>

using namespace NTL;

void keygen(ZZ& p_prime, ZZ& p, ZZ& n, ZZ& y, ZZ& D, const int l, const int k);

void precompute_pow2(ZZ& pow2k, const int k);

void encrypt(ZZ& c, const ZZ& m, const ZZ& n, const ZZ& y, const int k, const ZZ& pow2k);

void decrypt(ZZ& m, const ZZ& c, const ZZ& p, const ZZ& p_prime, const ZZ& D_, const int k, const ZZ& pow2k1);

void add_encrypted(ZZ& result, const ZZ& c1, const ZZ& c2, const ZZ& n);

void scalar_mult_encrypted(ZZ& result, const ZZ& c, const ZZ& scalar, const ZZ& n);

struct JLEncoding {
    Vec<ZZ> coeffs;
};

struct JLEncodingKeyPart {
    ZZ p_prime, p, n, g, D;
};

struct JLEncodingKey {
    Vec<JLEncodingKeyPart> keys;
    ZZ pow2k, pow2k1;
    long l, k, d;
};

JLEncodingKey gen_jl_encoding_key(long l, long k, long d);

JLEncoding encode(const ZZ_pE& m, const JLEncodingKey& key);

ZZ_pE decode(const JLEncoding& c, const JLEncodingKey& key);

void jle_add_assign(JLEncoding& a, const JLEncoding& b, const JLEncodingKey& key);

void jle_mult_assign(JLEncoding& a, const ZZ& scalar, const JLEncodingKey& key);