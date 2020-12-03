#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ_pE.h>
#include <NTL/ZZX.h>
#include <vector>

#include <gr.h>
#include <circuit.h>
#include <qrp.h>
#include <setup.h>
#include <rinocchio.h>

#include <fstream>
#include <sstream>
#include <string>

#include <assert.h>

using namespace std;
using namespace NTL;

void basicExample() {

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

Circuit circuitFromFile(string path) {
    Circuit c;
    ifstream File;
    File.open(path, ios::in);
    if (File) {
        File >> c;
    }
    return c;
}

void testMatrixMultCircuit(Circuit c) {
    ZZ modulus = ZZ(1) << 64;
    ZZ_p::init(modulus);

    ZZ_pX P = ZZ_pX();
    // P = x^4 + x + 1
    SetCoeff(P, 0);
    SetCoeff(P, 1);
    SetCoeff(P, 4);

    // P = x^32 + x^22 + x^2 + x^1 + 1
    // SetCoeff(P, 0);
    // SetCoeff(P, 1);
    // SetCoeff(P, 2);
    // SetCoeff(P, 22);
    // SetCoeff(P, 32);

    // instantiate GF(2^64, 4)
    ZZ_pE::init(P);
    cout << "modulus: " << P << "\n";

    const QRP qrp = getQRP(c);
    secretState state = setup(512, 64);
    CRS crs = getCRS(qrp, state);
    Vec<ZZ_p> input;
    input.SetLength(c.numberOfWires);
    input[0] = ZZ_p(1);
    for (long i = 1; i < c.numberOfWires; i++) {
        input[i] = ZZ_p(i);
    }

    Vec<ZZ_p> allWireValues = eval(c, input);
    cout << allWireValues << "all wires\n";
    Proof pi = prove(qrp, crs, allWireValues);

    Vec<ZZ_p> output;
    output.SetLength(c.numberOfOutputWires);
    for (long i = 0; i < c.numberOfOutputWires; i++) {
        output[i] = allWireValues[i + qrp.outOffset];
    }
    cout << output << "output\n";

    assert (verify(qrp, state, crs, pi, input, output) == 1);
}

int main() {
    basicExample();

    string path = "./out/matrix2.txt";
    path = "./out/formula.txt";
    const Circuit c = circuitFromFile(path);
    // printCircuit(c);
    testMatrixMultCircuit(c);
}