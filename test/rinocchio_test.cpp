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

#include <assert.h>

using namespace std;
using namespace NTL;

int main() {

    ZZ modulus = ZZ(1) << 64;
    ZZ_p::init(modulus);


    ZZ_pX P = ZZ_pX();
    // P = x^4 + x + 1
    // SetCoeff(P, 0);
    // SetCoeff(P, 1);
    // SetCoeff(P, 4);

    // P = x^32 + x^22 + x^2 + x^1 + 1
    SetCoeff(P, 0);
    SetCoeff(P, 1);
    SetCoeff(P, 22);
    SetCoeff(P, 32);

    // instantiate GF(2^64, 4)
    ZZ_pE::init(P);
    cout << "modulus: " << P << "\n";

    Circuit circuit;
    circuit.numberOfWires = 6;
    circuit.numberOfInputWires = 4;
    circuit.numberOfMidWires = 1;
    circuit.numberOfOutputWires = 1;
    circuit.numberOfMultiplicationGates = 2;
    Gate g4;
    g4.leftInputs = vector<long>();
    g4.leftInputs.push_back(2);
    g4.rightInputs = vector<long>();
    g4.rightInputs.push_back(3);
    Gate g5;
    g5.leftInputs = vector<long>();
    g5.leftInputs.push_back(0);
    g5.leftInputs.push_back(1);
    g5.rightInputs = vector<long>();
    g5.rightInputs.push_back(4);
    circuit.gates = vector<Gate>();
    circuit.gates.push_back(g4);
    circuit.gates.push_back(g5);

    QRP qrp = getQRP(circuit);
    secretState state = setup(512, 64);
    CRS crs = getCRS(qrp, state);
    Vec<ZZ_p> input;
    input.append(ZZ_p(3));
    input.append(ZZ_p(4));
    input.append(ZZ_p(2));
    input.append(ZZ_p(3));

    Vec<ZZ_p> allWireValues = eval(circuit, input);
    Proof pi = prove(qrp, crs, allWireValues);

    Vec<ZZ_p> output;
    output.append(ZZ_p(42));

    assert (verify(qrp, state, crs, pi, input, output) == 1);
}