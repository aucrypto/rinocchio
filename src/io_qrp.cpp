#include<io_qrp.h>
#include<gr.h>
#include<NTL/ZZ_pE.h>
#include<NTL/ZZX.h>

#include<time.h>
#include<fstream>
#include<sstream>
#include<algorithm>

using namespace NTL;
using namespace std;

ZZ_pEX ithDelta(const ZZ_pEX& targetPolynomial, const Vec<ZZ_pE>& exceptionalSet, long i) {
    ZZ_pE denominator;
    const ZZ_pEX x = conv<ZZ_pEX>("[[0] [1]]");
    const ZZ_pEX numerator = targetPolynomial / (x - exceptionalSet[i]);
    set(denominator);
    for (long j = 0; j < exceptionalSet.length(); j++)  {
        if (i == j) continue;
        denominator *= exceptionalSet[i] - exceptionalSet[j];
    }
    return numerator * getInverse(denominator);
    
}

IOQRP writeIOQRP(string path, const Circuit& circuit, long k, long minDegree) {
    IOQRP qrp;
    qrp.circuit = circuit;
    qrp.midOffset = circuit.numberOfInputWires;
    qrp.outOffset = circuit.numberOfWires - circuit.numberOfOutputWires;

    clock_t t = clock();
    const Vec<ZZ_pE> exceptionalSet = getExceptionalSubset(circuit.numberOfMultiplicationGates);
    BuildFromRoots(qrp.t, exceptionalSet);
    t = clock() - t;
    cout << "Target polynomial computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
    {
        ofstream File;
        File.open(path, ios::out);
        if (File) {
            File << qrp;
        } else {
            cout << "Could not write IO QRP to " << path << "\n";
        }
        File.close();
    }
    return qrp;
}

void writeYPolys(string path, const IOQRP& qrp) {
    clock_t t = clock();
    const Vec<ZZ_pE> exceptionalSet = getExceptionalSubset(qrp.circuit.numberOfMultiplicationGates);

    {   //Write y polynomial for all non input wires
        //y polynomial of input wires are implicitly zero
        ofstream File;
        File.open(path, ios::out);
        if (!File) {
            cout << "Could not write y " << path << "\n";
            return;
        }
        const ZZ_pEX x = conv<ZZ_pEX>("[[0] [1]]");
        ZZ_pE denominator;
        ZZ_pEX numerator;
        for (long i = 0; i < qrp.circuit.numberOfMultiplicationGates; i++) {
            div(numerator, qrp.t, (x - exceptionalSet[i]));
            set(denominator);
            for (long j = 0; j < exceptionalSet.length(); j++)  {
                if (i == j) continue;
                denominator *= exceptionalSet[i] - exceptionalSet[j];
            }
            mul(numerator, numerator, getInverse(denominator));
            File << numerator << endl;
        }
    }
    t = clock() - t;
    cout << "Computed and wrote deltas (ys) " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
    

}

void writeVandWPolys(string vPath, string wPath, string yPath, const IOQRP& qrp) {
    ZZ_pEX v_k, w_k, nextDelta;
    ofstream vFile, wFile;
    vFile.open(vPath, ios::out);
    wFile.open(wPath, ios::out);
    // qrp.t //figure out some nice way to get the max size of a line
    // long size = 10000000; //way to large
    // long size = 100000;
    // char line[size];
    string nextLine;
    for (long k = 0; k < qrp.circuit.numberOfWires; k++) {
        clear(v_k);
        clear(w_k);
        ifstream deltaFile(yPath);
        cout << k << "\n";
        for (long deltaIndex = 0; deltaIndex < qrp.circuit.numberOfMultiplicationGates; deltaIndex++) {
            // clock_t t = clock();
            // deltaFile >> nextDelta;
            // deltaFile.getline(line, size);
            // if (deltaFile.failbit) cout << "Read line failed\n";
            getline(deltaFile, nextLine);
            // t = clock() - t;
            // if ((deltaFile.rdstate() & std::ifstream::failbit ) != 0) {
            //     cout << "Read line failed. Buffer size: " << size << "\n";
            //     // todo iteratively attempt to fix buffer size?
            // }
            // cout << "time to read: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
            // t= clock();
            const vector<long> leftInputs = qrp.circuit.gates[deltaIndex].leftInputs;//todo use array?
            bool kthWireIsLeftInput = binary_search(leftInputs.begin(), leftInputs.end(), k);
            const vector<long> rightInputs = qrp.circuit.gates[deltaIndex].rightInputs;//todo use array?
            bool kthWireIsRightInput = binary_search(rightInputs.begin(), rightInputs.end(), k);
            if (kthWireIsLeftInput || kthWireIsRightInput) {
                istringstream ss;
                ss.str(nextLine);
                ss >> nextDelta;
                // deltaFile >> nextDelta;
                if (kthWireIsLeftInput) v_k += nextDelta;
                if (kthWireIsRightInput) w_k += nextDelta;
            } else {
                // deltaFile >> nextDelta;
                // deltaFile.ignore(numeric_limits<streamsize>::max(), '\n');
            }
        }
        vFile << v_k << endl;
        wFile << w_k << endl;
    }
}



ostream& operator<<(ostream& s, const IOQRP& qrp) {
    s << qrp.circuit;
    s << "\n";
    s << qrp.midOffset;
    s << "\n";
    s << qrp.outOffset;
    s << "\n";
    s << qrp.t;
    s << "\n";
    return s;
}
istream& operator>>(istream& s, IOQRP& qrp) {
    s >> qrp.circuit;
    s >> qrp.midOffset;
    s >> qrp.outOffset;
    s >> qrp.t;
    return s;
}
