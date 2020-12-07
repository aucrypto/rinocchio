#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
// #include <NTL/ZZ_limbs.h>
#include <assert.h>
#include <joye_libert.h>

#include <time.h>


using namespace std;
using namespace NTL;

void test() {
        ZZ modulus = ZZ(1) << 64;
    ZZ_p::init(modulus);
    
    // P = x^4 + x + 1
    ZZ_pX P = ZZ_pX();
    ZZ_p one = ZZ_p(1);
    SetCoeff(P, 0, one);
    SetCoeff(P, 1, one);
    SetCoeff(P, 4, one);

    // instantiate GF(2^64, 4)
    ZZ_pE::init(P);
    cout << "modulus: " << P << "\n";

    // TEST GENERATED IN SAGE:
    // Sage input:
    //     sage: R = Integers(2^64)
    //     sage: RX = PolynomialRing(R, 'x')
    //     sage: GR = QuotientRing(RX, RX.ideal(x^4 + x + 1))
    //     sage: GR
    //     Univariate Quotient Polynomial Ring in xbar over Ring of integers modulo 18446744073709551616 with modulus x^4 + x + 1
    //     sage: a = GR.random_element()
    //     sage: a
    //     2694002667186393808*xbar^3 + 550441077268878650*xbar^2 + 891207477139733263*xbar + 15182770761334213105
    //     sage: b = GR.random_element()
    //     sage: b
    //     8035703173926058224*xbar^3 + 3214717602632448704*xbar^2 + 17618849673318790798*xbar + 4043863594796898601
    //     sage: a+b
    //     10729705841112452032*xbar^3 + 3765158679901327354*xbar^2 + 63313076748972445*xbar + 779890282421560090
    //     sage: a*bw
    //     15044056529184173484*xbar^3 + 17074595192421832188*xbar^2 + 8353559141346242245*xbar + 15699992808455822249
    {
        ZZ_pE a = conv<ZZ_pE>("[15182770761334213105 891207477139733263 550441077268878650 2694002667186393808]");
        ZZ_pE b = conv<ZZ_pE>("[4043863594796898601 17618849673318790798 3214717602632448704 8035703173926058224]");
        ZZ_pE c_expected = conv<ZZ_pE>("[779890282421560090 63313076748972445 3765158679901327354 10729705841112452032]");
        ZZ_pE d_expected = conv<ZZ_pE>("[15699992808455822249 8353559141346242245 17074595192421832188 15044056529184173484]");
        ZZ_pE c_actual = a + b;
        ZZ_pE d_actual = a * b;

        assert (c_expected == c_actual);
        assert (d_expected == d_actual);
    }

    ZZ p_prime, p, n, y, D;
    const int k = 64; // message bit length
    const int l = 512; // modulus bit length

    keygen(p_prime, p, n, y, D, l, k);

    // cout << "p: " << p << "\n";
    // cout << "p_prime:" << p_prime << "\n";
    // cout << "n: " << n << "\n";
    // cout << "y: " << y << "\n";
    // cout << "D: " << D << "\n";

    ZZ pow2k, pow2k1;
    precompute_pow2(pow2k, k);
    precompute_pow2(pow2k1, k-1);

    ZZ m1 = RandomBits_ZZ(k-10); // so scalar multiplication doesn't exceed n
    ZZ m2 = RandomBits_ZZ(k);
    ZZ m1m2 = m1 + m2;
    ZZ s = ZZ(14);
    ZZ sm1 = s * m1;
    // cout << "m1:          " << m1 << "\n";
    // cout << "m2:          " << m2 << "\n";
    ZZ c1, c2, c1c2, sc1;
    encrypt(c1, m1, n, y, k, pow2k);
    encrypt(c2, m2, n, y, k, pow2k);
    add_encrypted(c1c2, c1, c2, n);
    scalar_mult_encrypted(sc1, c1, s, n);

    {
        // test homomorphic encryption
        ZZ m1_recovered, m2_recovered, m1m2_recovered, sm1_recovered;
        decrypt(m1_recovered, c1, p, p_prime, D, k, pow2k1);
        decrypt(m2_recovered, c2, p, p_prime, D, k, pow2k1);
        decrypt(m1m2_recovered, c1c2, p, p_prime, D, k, pow2k1);
        decrypt(sm1_recovered, sc1, p, p_prime, D, k, pow2k1);
        assert((m1 == m1_recovered) == 1);
        assert((m2 == m2_recovered) == 1);
        assert((m1m2 == m1m2_recovered) == 1);
        assert((sm1 == sm1_recovered) == 1);
    }

    {
        // test iterated addition
        ZZ m, m_sum, c, c_sum, m_sum_decrypted;
        m_sum = RandomBits_ZZ(k);
        encrypt(c_sum, m_sum, n, y, k, pow2k);
        for (int i = 0; i < 200; i++) {
            m = RandomBits_ZZ(k);
            encrypt(c, m, n, y, k, pow2k);
            m_sum = (m_sum + m) % modulus;
            add_encrypted(c_sum, c_sum, c, n);
        }
        decrypt(m_sum_decrypted, c_sum, p, p_prime, D, k, pow2k1);
        // cout << "c_sum          : " << c_sum << "\n";
        // cout << "m_sum          : " << m_sum << "\n";
        // cout << "m_sum_decrypted: " << m_sum_decrypted << "\n";
        assert((m_sum == m_sum_decrypted) == 1);
    }

    {
        // test iterated scaling
        ZZ scalar, m_scaled, c, c_scaled, m_scaled_decrypted;
        m_scaled = RandomBits_ZZ(k);
        encrypt(c_scaled, m_scaled, n, y, k, pow2k);
        for (int i = 0; i < 200; i++) {
            scalar = RandomBits_ZZ(k);
            m_scaled = (m_scaled * scalar) % modulus;
            scalar_mult_encrypted(c_scaled, c_scaled, scalar, n);
        }
        decrypt(m_scaled_decrypted, c_scaled, p, p_prime, D, k, pow2k1);
        // cout << "c_sum          : " << c_sum << "\n";
        // cout << "m_sum          : " << m_sum << "\n";
        // cout << "m_sum_decrypted: " << m_sum_decrypted << "\n";
        assert((m_scaled == m_scaled_decrypted) == 1);
    }

    {
        ZZ_pE m1 = random_ZZ_pE();
        ZZ_pE m2 = random_ZZ_pE();
        ZZ_pE m3 = random_ZZ_pE();
        JLEncodingKey key = gen_jl_encoding_key(l, k);
        JLEncoding encoded1 = encode(m1, key);
        JLEncoding encoded2 = encode(m2, key);
        JLEncoding encoded3 = encode(m3, key);
        
        ZZ_pE decoded = decode(encoded1, key);
        // cout << "m1......: " << m1 << "\n";
        // cout << "decoded: " << decoded << "\n";
        assert (m1 == decoded);

        JLEncoding encoded0;
        jle_add_assign(encoded0, encoded1, key);
        jle_add_assign(encoded1, encoded2, key);
        jle_add_assign(encoded1, encoded3, key);
        decoded = decode(encoded0, key);
        assert (decoded == (m1));
        decoded = decode(encoded1, key);
        assert (decoded == (m1 + m2 + m3));
        
        ZZ scalar = ZZ(1231314);
        jle_scalar_mult_assign(encoded2, scalar, key);
        decoded = decode(encoded2, key);
        assert (decoded == (conv<ZZ_pE>(scalar) * m2));

    }

    
    {
        ZZ_pE m1 = random_ZZ_pE();
        ZZ_pE m2 = random_ZZ_pE();
        JLEncodingKey key = gen_jl_encoding_key(l, k);
        JLEncoding encoded1 = encode(m1, key);
        JLEncoding encoded2 = encode(m2, key);
        
        ZZ_pE decoded = decode(encoded1, key);
        // cout << "m1......: " << m1 << "\n";
        // cout << "decoded: " << decoded << "\n";
        assert (m1 == decoded);

        JLEncoding encoded12prod = jle_mult(encoded1, to_vec_ZZ(rep(m2).rep), key);
        decoded = decode(encoded12prod, key);
        assert (decoded == (m1 * m2));

        JLEncoding encoded12prod_ = PlainMulEncryption(encoded1, to_vec_ZZ(rep(m2).rep), key.n);
        ZZ_pE decoded2 = decode(encoded12prod_, key);
        assert (decoded2 == (m1*m2));
    }

    cout << "\nPassed all tests!\n";
}

int main() {
    test();

    long iterations = 100000;
    cout << "Testing with " << iterations << " iterations\n";

    clock_t t = clock();
    ZZ modulus = ZZ(1) << 64;


    ZZ p_prime, p, n, y, D;
    const int k = 64; // message bit length
    const int l = 512; // modulus bit length

    keygen(p_prime, p, n, y, D, l, k);

    ZZ pow2k, pow2k1;
    precompute_pow2(pow2k, k);
    precompute_pow2(pow2k1, k-1);

    Vec<ZZ> randomMods, encryptions;
    randomMods.SetLength(iterations);
    encryptions.SetLength(iterations);

    t = clock();
    for (int i = 0; i < iterations; i++) {
        randomMods[i] = RandomBits_ZZ(k);
    }    
    t = clock() - t;
    cout << "Draw random numbers mod 2^k: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    t = clock();
    ZZ m, c;
    for (int i = 0; i < iterations; i++) {
        encrypt(c, randomMods[i], n, y, k, pow2k);
        encryptions[i] = c;
    }    
    t = clock() - t;
    cout << "Encryption: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    {
        // test iterated addition
        ZZ m, m_sum, c, c_sum, m_sum_decrypted;
        m_sum = RandomBits_ZZ(k);
        t = clock();
        encrypt(c_sum, m_sum, n, y, k, pow2k);
        for (int i = 0; i < iterations; i++) {
            m = randomMods[i];
            m_sum = m_sum + m;
            add_encrypted(c_sum, c_sum, encryptions[i], n);
        }
        decrypt(m_sum_decrypted, c_sum, p, p_prime, D, k, pow2k1);
        t = clock() - t;
        cout << "Addition: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
        assert((m_sum % modulus) == m_sum_decrypted);
    }

    {
        // test iterated scaling
        ZZ scalar, m_scaled, c, c_scaled, m_scaled_decrypted;
        m_scaled = RandomBits_ZZ(k);
        encrypt(c_scaled, m_scaled, n, y, k, pow2k);
        t = clock();
        for (int i = 0; i < iterations; i++) {
            scalar = randomMods[i];
            MulMod(m_scaled, m_scaled, scalar, modulus);
            scalar_mult_encrypted(c_scaled, c_scaled, scalar, n);
        }
        decrypt(m_scaled_decrypted, c_scaled, p, p_prime, D, k, pow2k1);
        t = clock() - t;
        cout << "Scalar multiplication: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
        assert(m_scaled == m_scaled_decrypted);
    }

    {
        // return 0;
        // test iterated decryption
        Vec<ZZ> decryptions;
        decryptions.SetLength(iterations);
        t = clock();
        for (int i = 0; i < iterations; i++) {
            decrypt(decryptions[i], encryptions[i], p, p_prime, D, k, pow2k1);
        }
        t = clock() - t;
        cout << "Decryption: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
        assert(randomMods == decryptions);
    }
}