#include<qrp.h>
#include<gr.h>
#include<NTL/ZZ_pE.h>
#include<NTL/ZZX.h>

using namespace NTL;
using namespace std;

Vec<ZZ_pEX> getInterpolationDeltas(long numberOfMultiplicationGates, ZZ_pEX& targetPolynomial) {
    Vec<ZZ_pEX> deltas;
    deltas.SetLength(numberOfMultiplicationGates);

    ZZ_pEX x;
    SetX(x);
    set(targetPolynomial);
    cout << "Computing subset of exceptional set\n";
    const Vec<ZZ_pE> exceptionalSubSet = getExceptionalSubset(numberOfMultiplicationGates);
    cout << "Subset of exceptional set computed\n";
    const Vec<ZZ_pEX> termsOfT = getTargetPolynomialTerms(numberOfMultiplicationGates);
    cout << "terms of t computed\n"; //We get here very fast
    // Vec<ZZ_pE> exceptionalSubSet;
    // exceptionalSubSet.SetLength(numberOfMultiplicationGates);
    for (long i = 0; i < numberOfMultiplicationGates; i++)  {
        // cout << "Compute exceptional elements and target: " << i << endl;
        // const ZZ_pE exceptionalElement = exceptionalSubSet[i];
        // targetPolynomial *= x - exceptionalElement;
        targetPolynomial *= termsOfT[i];//todo I think this was slower
    }
    cout << "Target polynomial computed\n";//slow, would using a c style array be any faster?

    ZZ_pE denominator;
    for (long i = 0; i < numberOfMultiplicationGates; i++)  {
        // cout << "Compute deltas: " << i << endl;
        const ZZ_pEX numerator = targetPolynomial / (x - exceptionalSubSet[i]);
        set(denominator);
        for (long j = 0; j < numberOfMultiplicationGates; j++)  {
            if (i == j) continue;

            denominator *= exceptionalSubSet[i] - exceptionalSubSet[j];
        }
        deltas[i] = numerator * getInverse(denominator);
    }
    cout << "Deltas computed\n";//slow

    return deltas;
}

QRP getQRP(const Circuit& circuit) {
    //Pick distinct elements of exceptional set for each gate:
    ZZ_pEX targetPolynomial;
    Vec<ZZ_pEX> deltas = getInterpolationDeltas(circuit.numberOfMultiplicationGates, targetPolynomial);

    Vec<ZZ_pEX> V, W, Y;
    V.SetLength(circuit.numberOfWires);
    W.SetLength(circuit.numberOfWires);
    Y.SetLength(circuit.numberOfWires);

    ZZ_pEX x;
    SetX(x);
    
    for (int k = 0; k < circuit.numberOfMultiplicationGates; k++) {

        int nextLeftInput = 0;
        int nextRightInput = 0;
        vector<long> leftInputs = circuit.gates[k].leftInputs;
        vector<long> rightInputs = circuit.gates[k].rightInputs;
        for (int j = 0; j < circuit.numberOfWires; j++) {
            if (j == k + circuit.numberOfInputWires) {
                Y[j] += deltas[k];
                // cout << j << " is output of " << k + circuit.numberOfInputWires << endl;
            }
            
            if (nextLeftInput < leftInputs.size() && j == leftInputs[nextLeftInput]) {
                // cout << j << " is left input of " << k + circuit.numberOfInputWires << endl;
                V[j] += deltas[k];
                nextLeftInput++;
            } 
            
            if (nextRightInput < rightInputs.size() && j == rightInputs[nextRightInput]) {
                // cout << j << " is right input of " << k + circuit.numberOfInputWires << endl;
                W[j] += deltas[k];
                nextRightInput++;
            }
        }
    }
    cout << "V, W and Y computed\n";
    
    QRP qrp;
    qrp.circuit = circuit;
    qrp.midOffset = circuit.numberOfInputWires;
    qrp.outOffset = circuit.numberOfWires - circuit.numberOfOutputWires;
    qrp.t = targetPolynomial;
    qrp.V = V;
    qrp.W = W;
    qrp.Y = Y;
    cout << "qrp function exit\n";
    return qrp;
}


ostream& operator<<(ostream& s, const QRP& qrp) {
    s << qrp.circuit;
    s << "\n";
    s << qrp.midOffset;
    s << "\n";
    s << qrp.outOffset;
    s << "\n";
    s << qrp.t;
    s << "\n";
    s << qrp.V;
    s << "\n";
    s << qrp.W;
    s << "\n";
    s << qrp.Y;
    s << "\n";

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
