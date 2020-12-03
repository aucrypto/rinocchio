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
    Vec<ZZ_pE> exceptionalSubSet;
    exceptionalSubSet.SetLength(numberOfMultiplicationGates);
    for (long i = 0; i < numberOfMultiplicationGates; i++)  {
        cout << "Compute exceptional elements and target: " << i << endl;
        const ZZ_pE exceptionalElement = indexedElementInExceptionalSet(i);
        exceptionalSubSet[i] = exceptionalElement;
        targetPolynomial *= x - exceptionalElement;
    }

    ZZ_pE denominator;
    for (long i = 0; i < numberOfMultiplicationGates; i++)  {
        cout << "Compute deltas: " << i << endl;
        const ZZ_pEX numerator = targetPolynomial / (x - exceptionalSubSet[i]);
        set(denominator);
        for (long j = 0; j < numberOfMultiplicationGates; j++)  {
            if (i == j) continue;

            denominator *= exceptionalSubSet[i] - exceptionalSubSet[j];
        }
        deltas[i] = numerator * getInverse(denominator);
    }

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
        cout << "Compute " << k << endl;

        Gate kthMultGate =  circuit.gates[k];

        int nextLeftInput = 0;
        int nextRightInput = 0;
        vector<long> leftInputs = kthMultGate.leftInputs;
        vector<long> rightInputs = kthMultGate.rightInputs;
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
    
    QRP qrp;
    qrp.circuit = circuit;
    qrp.midOffset = circuit.numberOfInputWires;
    qrp.outOffset = circuit.numberOfWires - circuit.numberOfOutputWires;
    qrp.t  = targetPolynomial;
    qrp.V = V;
    qrp.W = W;
    qrp.Y = Y;
    return qrp;
}
