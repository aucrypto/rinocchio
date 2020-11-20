#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>

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
    //     sage: a*b
    //     15044056529184173484*xbar^3 + 17074595192421832188*xbar^2 + 8353559141346242245*xbar + 15699992808455822249

    cout << P << "\n";
}