#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
// #include <NTL/ZZ_limbs.h>
#include <assert.h>
#include "../include/joy_libert.h"


using namespace std;
using namespace NTL;

int main() {
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

    cout << "p: " << p << "\n";
    cout << "p_prime:" << p_prime << "\n";
    cout << "n: " << n << "\n";
    cout << "y: " << y << "\n";
    cout << "D: " << D << "\n";

    ZZ pow2k, pow2k1;
    precompute_pow2(pow2k, k);
    precompute_pow2(pow2k1, k-1);

    ZZ m = RandomBits_ZZ(k);
    cout << "m:           " << m << "\n";
    ZZ c;
    encrypt(c, m, n, y, k, pow2k);
    ZZ m_recovered;
    decrypt(m_recovered, c, p, p_prime, D, k, pow2k1);

    cout << "m_recovered: " << m_recovered << "\n";

    assert((m == m_recovered) == 1);

    cout << "\nPassed all tests!\n";
}