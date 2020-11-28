#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ_pE.h>
#include <NTL/ZZX.h>
#include <vector>

#include <gr.h>
#include <qrp.h>
#include <setup.h>
#include <rinocchio.h>

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

    Circuit circuit;
    circuit.numberOfWires = 6;
    circuit.numberOfInputWires = 4;
    circuit.numberOfMidWires = 1;
    circuit.numberOfOutputWires = 1;
    circuit.numberOfMultiplicationGates = 2;

    QRP qrp = getQRP(circuit);
    secretState state = setup();
    CRS crs = getCRS(qrp, state);
    Vec<ZZ_p> input;
    input.append(ZZ_p(3));
    input.append(ZZ_p(4));
    input.append(ZZ_p(2));
    input.append(ZZ_p(3));

    Vec<ZZ_p> allWireValues = input;
    //Compute values for all multiplication gates:
    ZZ_p c_5, c_6;//todo generalize using circuit representation
    
    c_5 = input(3) * input(4);
    allWireValues.append(c_5);
    cout << c_5 << " should be 6\n";

    c_6 = c_5 * (input(1) + input(2));
    allWireValues.append(c_6);
    cout << c_6 << " should be 42\n";
    Proof pi = prove(qrp, crs, allWireValues);

    Vec<ZZ_p> output;
    output.append(ZZ_p(42));

    cout << verify(qrp, state, crs, pi, input, output) << "\n";
}