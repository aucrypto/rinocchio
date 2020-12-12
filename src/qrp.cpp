#include<qrp.h>
#include<gr.h>
#include<NTL/ZZ_pE.h>
#include<NTL/ZZX.h>

#include<time.h>

using namespace NTL;
using namespace std;

ZZ_pEX kthDelta(const ZZ_pEX& targetPolynomial, const Vec<ZZ_pE>& exceptionalSet, long i) {
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

QRP getQRP(const Circuit& circuit, long k, long minDegree) {
    QRP qrp;
    qrp.circuit = circuit;
    qrp.midOffset = circuit.numberOfInputWires;
    qrp.outOffset = circuit.numberOfWires - circuit.numberOfOutputWires;
    // ZZ mod = ZZ(1) << k;
    // ZZ_p::init(mod);
    // ZZ_pX polyMod = primitiveIrredPoly(minDegree); //todo check if minDegree is >= log_2(multGates)
    // ZZ_pE::init(polyMod);
    // todo fix automatic minimal modulus
    
    //Pick distinct elements of exceptional set for each gate:
    cout << "Computing subset of exceptional set...\n";
    clock_t t = clock();
    const Vec<ZZ_pE> exceptionalSet = getExceptionalSubset(circuit.numberOfMultiplicationGates);
    t = clock() - t;
    cout << "Subset of exceptional set computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    t = clock();
    BuildFromRoots(qrp.t, exceptionalSet);
    t = clock() - t;
    cout << "Target polynomial computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    qrp.V.SetLength(circuit.numberOfWires);
    qrp.W.SetLength(circuit.numberOfWires);
    qrp.Y.SetLength(circuit.numberOfWires);
    

    const ZZ_pEX x = conv<ZZ_pEX>("[[0] [1]]");
    ZZ_pE denominator;
    ZZ_pEX numerator;
    t = clock();
    for (int k = 0; k < circuit.numberOfMultiplicationGates; k++) {
        div(numerator, qrp.t, (x - exceptionalSet[k]));
        set(denominator);
        for (long j = 0; j < exceptionalSet.length(); j++)  {
            if (k == j) continue;
            denominator *= exceptionalSet[k] - exceptionalSet[j];
        }
        mul(numerator, numerator, getInverse(denominator));
        // ZZ_pEX delta = kthDelta(qrp.t, exceptionalSubSet, k);
        qrp.Y[k + circuit.numberOfInputWires] = numerator;
        int nextLeftInput = 0;
        int nextRightInput = 0;
        vector<long> leftInputs = circuit.gates[k].leftInputs;
        vector<long> rightInputs = circuit.gates[k].rightInputs;
        for (int j = 0; j < circuit.numberOfWires; j++) {
            // if (j == k + circuit.numberOfInputWires) {
            //     qrp.Y[j] += delta;
            //     // cout << j << " is output of " << k + circuit.numberOfInputWires << endl;
            // }
            
            if (nextLeftInput < leftInputs.size() && j == leftInputs[nextLeftInput]) {
                // cout << j << " is left input of " << k + circuit.numberOfInputWires << endl;
                qrp.V[j] += numerator;
                nextLeftInput++;
            } 
            
            if (nextRightInput < rightInputs.size() && j == rightInputs[nextRightInput]) {
                // cout << j << " is right input of " << k + circuit.numberOfInputWires << endl;
                qrp.W[j] += numerator;
                nextRightInput++;
            }
        }
    }
    t = clock() - t;
    cout << "V, W and Y computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
    return qrp;
}

void getLargeQRP(const Circuit& circuit, long k, long minDegree) {
    QRP qrp;
    qrp.circuit = circuit;
    qrp.midOffset = circuit.numberOfInputWires;
    qrp.outOffset = circuit.numberOfWires - circuit.numberOfOutputWires;
    // ZZ mod = ZZ(1) << k;
    // ZZ_p::init(mod);
    // ZZ_pX polyMod = primitiveIrredPoly(minDegree); //todo check if minDegree is >= log_2(multGates)
    // ZZ_pE::init(polyMod);
    // todo fix automatic minimal modulus
    
    //Pick distinct elements of exceptional set for each gate:
    cout << "Computing subset of exceptional set...\n";
    clock_t t = clock();
    const Vec<ZZ_pE> exceptionalSet = getExceptionalSubset(circuit.numberOfMultiplicationGates);
    t = clock() - t;
    cout << "Subset of exceptional set computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    t = clock();
    BuildFromRoots(qrp.t, exceptionalSet);
    t = clock() - t;
    cout << "Target polynomial computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    qrp.V.SetLength(circuit.numberOfWires);
    qrp.W.SetLength(circuit.numberOfWires);
    qrp.Y.SetLength(circuit.numberOfWires);
    

    const ZZ_pEX x = conv<ZZ_pEX>("[[0] [1]]");
    ZZ_pE denominator;
    ZZ_pEX numerator;
    t = clock();
    for (int k = 0; k < circuit.numberOfMultiplicationGates; k++) {
        div(numerator, qrp.t, (x - exceptionalSet[k]));
        set(denominator);
        for (long j = 0; j < exceptionalSet.length(); j++)  {
            if (k == j) continue;
            denominator *= exceptionalSet[k] - exceptionalSet[j];
        }
        mul(numerator, numerator, getInverse(denominator));
        // ZZ_pEX delta = kthDelta(qrp.t, exceptionalSubSet, k);
        qrp.Y[k + circuit.numberOfInputWires] = numerator;
        int nextLeftInput = 0;
        int nextRightInput = 0;
        vector<long> leftInputs = circuit.gates[k].leftInputs;
        vector<long> rightInputs = circuit.gates[k].rightInputs;
        for (int j = 0; j < circuit.numberOfWires; j++) {
            // if (j == k + circuit.numberOfInputWires) {
            //     qrp.Y[j] += delta;
            //     // cout << j << " is output of " << k + circuit.numberOfInputWires << endl;
            // }
            
            if (nextLeftInput < leftInputs.size() && j == leftInputs[nextLeftInput]) {
                // cout << j << " is left input of " << k + circuit.numberOfInputWires << endl;
                qrp.V[j] += numerator;
                nextLeftInput++;
            } 
            
            if (nextRightInput < rightInputs.size() && j == rightInputs[nextRightInput]) {
                // cout << j << " is right input of " << k + circuit.numberOfInputWires << endl;
                qrp.W[j] += numerator;
                nextRightInput++;
            }
        }
    }
    t = clock() - t;
    cout << "V, W and Y computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
}


ostream& operator<<(ostream& s, const QRP& qrp) {
    clock_t t = clock();
    s << qrp.circuit;
    s << "\n";
    s << qrp.midOffset;
    s << "\n";
    s << qrp.outOffset;
    s << "\n";
    t = clock() - t;
    cout << "Wrote QRP circuit: " << ((double) t) / CLOCKS_PER_SEC << " seconds.\n";
    t = clock();
    s << qrp.t;
    s << "\n";
    t = clock() - t;
    cout << "Wrote QRP target: " << ((double) t) / CLOCKS_PER_SEC << " seconds.\n";
    t = clock();
    s << qrp.V;
    s << "\n";
    t = clock() - t;
    cout << "Wrote QRP V: " << ((double) t) / CLOCKS_PER_SEC << " seconds.\n";
    t = clock();
    s << qrp.W;
    s << "\n";
    t = clock() - t;
    cout << "Wrote QRP W: " << ((double) t) / CLOCKS_PER_SEC << " seconds.\n";
    t = clock();
    s << qrp.Y;
    s << "\n";
    t = clock() - t;
    cout << "Wrote QRP Y: " << ((double) t) / CLOCKS_PER_SEC << " seconds.\n";

    return s;
}
istream& operator>>(istream& s, QRP& qrp) {
    s >> qrp.circuit;
    // long next;
    // if (!(s >> next)) {
    //     cout << "not\n";
    // }
    //todo handle errors
    s >> qrp.midOffset;
    s >> qrp.outOffset;
    s >> qrp.t;
    s >> qrp.V;
    s >> qrp.W;
    s >> qrp.Y;

    return s;
}
