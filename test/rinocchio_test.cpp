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
#include <io_qrp.h>
#include <setup.h>
#include <rinocchio.h>

#include <fstream>
#include <sstream>
#include <string>

#include <assert.h>

#include <time.h>

using namespace std;
using namespace NTL;


void init(long k, long degree) {

    ZZ modulus = ZZ(1) << k;
    ZZ_p::init(modulus);

    ZZ_pX P = primitiveIrredPoly(degree);
    ZZ_pE::init(P);
    cout << "modulus: " << P << "\n";
}

void testMatrixMultCircuit(const QRP& qrp);

void basicExample() {

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
    // SetCoeff(P, 22);
    // SetCoeff(P, 32);

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

    QRP qrp = getQRP(circuit, 64, 4);
    SecretState state = setup(qrp, 512, 64);
    CRS crs = getCRS(qrp, state);
    Vec<ZZ_p> input;
    input.append(ZZ_p(3));
    input.append(ZZ_p(4));
    input.append(ZZ_p(2));
    input.append(ZZ_p(3));

    Vec<ZZ_p> allWireValues = eval(circuit, input);
    Proof pi = prove(qrp, crs, allWireValues);

    Vec<ZZ_p> output;
    output.append(allWireValues[5]);

    assert (verify(state, crs, pi, input, output) == 1);
}

Circuit circuitFromFile(string path) {
    Circuit c;
    ifstream File;
    File.open(path, ios::in);
    if (File) {
        File >> c;
        File.close();
    } else {
        cout << "not file\n";
    }
    return c;
}

void writeQRP(const QRP& qrp, string path) {
    ofstream File;
    File.open(path, ios::out);
    if (File) {
        File << qrp;
    } else {
        cout << "not file\n";
    }
}

QRP qrpFromFile(string path, bool& success) {
    QRP qrp;
    ifstream File;
    File.open(path, ios::in);
    if (File) {
        success = true;
        File >> qrp;
        File.close();
    } else {
        success = false;
    }
    return qrp;
}

void writeSecretState(const SecretState& ss, string path) {
    ofstream File;
    File.open(path, ios::out);
    if (File) {
        File << ss;
    } else {
        cout << "not file\n";
    }
}

SecretState readSecretState(string path, bool& success) {
    SecretState secretState;
    ifstream File;
    File.open(path, ios::in);
    if (File) {
        success = true;
        File >> secretState;
        File.close();
    } else {
        success = false;
    }
    return secretState;
}

void writeCRS(const CRS& crs, string path) {
    ofstream File;
    File.open(path, ios::out);
    if (File) {
        File << crs;
    } else {
        cout << "not file\n";
    }
}

CRS readCRS(string path, bool& success) {
    CRS crs;
    ifstream File;
    File.open(path, ios::in);
    if (File) {
        success = true;
        File >> crs;
        File.close();
    } else {
        success = false;
    }
    return crs;
}

QRP computeOrReadQRP(string qrpPath, string circuitPath) {
    bool readSuccess;
    QRP qrp = qrpFromFile(qrpPath, readSuccess);
    if (!readSuccess) {
        cout << "QRP file does not exist.\n";
        Circuit c = circuitFromFile(circuitPath); //todo handle error
        qrp = getQRP(c, 64, 10);
        writeQRP(qrp, qrpPath);
        cout << "QRP was computed and written to file.\n";
        return qrp;
    }
    if (qrp.circuit.numberOfWires == 0
        || qrp.circuit.numberOfMultiplicationGates != qrp.circuit.gates.size()
        || qrp.circuit.numberOfWires != qrp.V.length()
        || qrp.circuit.numberOfWires != qrp.W.length()
        || qrp.circuit.numberOfWires != qrp.Y.length()) {
            cout << "Read invalid QRP\n"; //todo maybe ask before overwriting
            Circuit c = circuitFromFile(circuitPath); //todo handle error
            qrp = getQRP(c, 64, 10);
            writeQRP(qrp, qrpPath);
            cout << "QRP was computed and written to file.\n";
            return qrp;
    }
    cout << "QRP was read from file.\n";
    return qrp;
}

void testMatrixMultCircuit(const QRP& qrp, const SecretState& state, const CRS& crs) {
    Vec<ZZ_p> input;
    input.SetLength(qrp.circuit.numberOfInputWires);
    input[0] = ZZ_p(1);
    for (long i = 1; i < qrp.circuit.numberOfInputWires; i++) {
        input[i] = to_ZZ_p(RandomBits_ZZ(64));
    }

    clock_t t = clock();
    Vec<ZZ_p> allWireValues = eval(qrp.circuit, input);
    t = clock() - t;
    cout << "Circuit evaluated: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    Vec<ZZ_p> output;
    output.SetLength(qrp.circuit.numberOfOutputWires);
    for (long i = 0; i < qrp.circuit.numberOfOutputWires; i++) {
        output[i] = allWireValues[i + qrp.outOffset];
    }

    t = clock();
    Proof pi = prove(qrp, crs, allWireValues);
    t = clock() - t;
    cout << "Proof done: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    t = clock();
    assert (verify(state, crs, pi, input, output) == 1);
    t = clock() - t;
    cout << "Verify done: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
}

void testFile(string testName, long k) {
    string outDir = "./out/setups/" + testName + "/";
    string path = outDir + testName;
    string circuitPath = path + "_circuit.txt";
    string secretPath = path + "_secret.txt";
    string crsPath = path + "_crs.txt";
    string qrpPath = path + "_qrp.txt";
    string ioqrpPath = path  + "_io_qrp.txt";
    const Circuit c = circuitFromFile(circuitPath);

    long extensionDegree = 2;
    if (c.numberOfMultiplicationGates > 4) extensionDegree = 3;
    if (c.numberOfMultiplicationGates > 8) extensionDegree = 4;
    if (c.numberOfMultiplicationGates > 16) extensionDegree = 5;
    if (c.numberOfMultiplicationGates > 32) extensionDegree = 6;
    if (c.numberOfMultiplicationGates > 64) extensionDegree = 7;
    if (c.numberOfMultiplicationGates > 128) extensionDegree = 8;
    if (c.numberOfMultiplicationGates > 256) extensionDegree = 9;
    if (c.numberOfMultiplicationGates > 512) extensionDegree = 10;
    if (c.numberOfMultiplicationGates > 1024) extensionDegree = 11;
    if (c.numberOfMultiplicationGates > 2048) extensionDegree = 12;
    if (c.numberOfMultiplicationGates > 4096) extensionDegree = 13;
    if (c.numberOfMultiplicationGates > 8192) extensionDegree = 14;
    if (c.numberOfMultiplicationGates > 16384) extensionDegree = 15;
    if (c.numberOfMultiplicationGates > 32768) extensionDegree = 16;
    init(k, extensionDegree);

    

    // return;
    const QRP qrp = computeOrReadQRP(qrpPath, circuitPath);


    clock_t t = clock();
    bool success;
    SecretState state = readSecretState(secretPath, success);
    if (!success) {
        state = setup(qrp, 512, 64);
        writeSecretState(state, secretPath);
    }
    t = clock() - t;
    cout << "Secret done: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
    
    t = clock();
    CRS crs = readCRS(crsPath, success);
    if (!success) {
        crs = getCRS(qrp, state);
        writeCRS(crs, crsPath);
    }
    t = clock() - t;
    cout << "CRS done: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
    
    testMatrixMultCircuit(qrp, state, crs);

}

void testnxnxnMatrixMult(int size) {
    string testName = string("n=") + to_string(size) +  "_m=" + to_string(size) + "_k=" + to_string(size);
    cout << "Testing: " << testName << "\n";
    clock_t t = clock();
    testFile(testName, 64);
    t = clock() - t;
    cout << "Tested: " << testName << "\n";
    cout << "Spent " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
}

/**
 * Instantiates the protocol with GR(2^k, d) to achieve soundness error 2^-d
 * Precondtion: d is one of 40, 60, or 80
 **/
void testNaiveSoundness(int size, long k, long d) {
    string testName = string("n=") + to_string(size) +  "_m=" + to_string(size) + "_k=" + to_string(size);
    cout << "Testing: " << testName << ", soundness error 2^-" << d << "\n";
    init(k, d);

    string outDir = "./out/setups/" + testName + "/";
    string path = outDir + testName;
    string circuitPath = path + "_circuit.txt";
    const Circuit c = circuitFromFile(circuitPath);

    string qrpPath = path + "_d=" + to_string(d) + "_qrp.txt";
    const QRP qrp = computeOrReadQRP(qrpPath, circuitPath);

    string secretPath = path + "_d=" + to_string(d) + "_secret.txt";
    string crsPath = path + "_d=" + to_string(d) + "_crs.txt";

    clock_t t = clock();
    bool success;
    SecretState state = readSecretState(secretPath, success);
    if (!success) {
        cout << "Computing secret state\n";
        t = clock();
        state = setup(qrp, 512, 64);
        t = clock() - t;
        cout << "Secret done: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
        writeSecretState(state, secretPath);
    }
    
    CRS crs = readCRS(crsPath, success);
    if (!success) {
        cout << "Computing CRS\n";
        t = clock();
        crs = getCRS(qrp, state);
        t = clock() - t;
        cout << "CRS done: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
        writeCRS(crs, crsPath);
    }
    
    testMatrixMultCircuit(qrp, state, crs);

}

int main() {
    // testNaiveSoundness(10, 64, 40);
    // return 0;
    cout << "##################################################" << endl;
    cout << "--------------------------40----------------------" << endl;
    cout << "##################################################" << endl;
    for (int i = 2; i <= 10; i++) {

        testNaiveSoundness(i, 64, 40);
    }

    cout << "##################################################" << endl;
    cout << "--------------------------60----------------------" << endl;
    cout << "##################################################" << endl;

    for (int i = 2; i <= 10; i++) {

        testNaiveSoundness(i, 64, 60);
    }

    cout << "##################################################" << endl;
    cout << "--------------------------80----------------------" << endl;
    cout << "##################################################" << endl;

    for (int i = 2; i <= 10; i++) {

        testNaiveSoundness(i, 64, 80);
    }
    // testnxnxnMatrixMult(10);
    // return 0;
    // for (int i = 2; i <=25; i++) {
    //     testnxnxnMatrixMult(i);
    // }
}