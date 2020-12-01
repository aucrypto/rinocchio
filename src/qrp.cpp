#include<qrp.h>
#include<gr.h>
#include<NTL/ZZ_pE.h>
#include<NTL/ZZX.h>

using namespace NTL;
using namespace std;

Vec<ZZ_pEX> getInterpolationDeltas(long numberOfMultiplicationGates) {
    Vec<ZZ_pEX> deltas;
    deltas.SetLength(numberOfMultiplicationGates);

    //todo we could generate all the needed elements of the exceptional set only once, and more efficiently

    ZZ_pEX x, numerator;
    ZZ_pE denominator;
    SetX(x);
    for (long i = 0; i < numberOfMultiplicationGates; i++)  {
        set(numerator);
        set(denominator);
        for (long j = 0; j < numberOfMultiplicationGates; j++)  {
            if (i == j) continue;
            numerator *= x - indexedElementInExceptionalSet(j);
            denominator *= indexedElementInExceptionalSet(i) - indexedElementInExceptionalSet(j);
        }
        deltas[i] = numerator * getInverse(denominator);
    }

    return deltas;
}

QRP getQRP(Circuit circuit) {
    //Pick distinct elements of exceptional set for each gate:
    Vec<ZZ_pEX> deltas = getInterpolationDeltas(circuit.numberOfMultiplicationGates);

    Vec<ZZ_pEX> V, W, Y;
    V.SetLength(circuit.numberOfWires);
    W.SetLength(circuit.numberOfWires);
    Y.SetLength(circuit.numberOfWires);

    // Compute t(x) = (x - g_1) * (x - g_2) ...
    ZZ_pEX t = conv<ZZ_pEX>(1);
    ZZ_pEX x;
    SetX(x);
    
    for (int k = 0; k < circuit.numberOfMultiplicationGates; k++) {

        //target polynomial:
        t *= x - indexedElementInExceptionalSet(k); // this is used when generating deltas as well

        Gate kthMultGate =  circuit.gates[k];

        int nextLeftInput = 0;
        int nextRightInput = 0;
        vector<long> leftInputs = kthMultGate.leftInputs;
        vector<long> rightInputs = kthMultGate.rightInputs;
        for (int j = 0; j < circuit.numberOfWires; j++) {
            if (j == k + circuit.numberOfInputWires) {
                Y[j] += deltas[k];
            }
            
            if (nextLeftInput < leftInputs.size() && j == leftInputs[nextLeftInput]) {
                V[j] += deltas[k];
                nextLeftInput++;
            } 
            
            if (nextRightInput < rightInputs.size() && j == rightInputs[nextRightInput]) {
                W[j] += deltas[k];
                nextRightInput++;
            }
        }
    }
    
    QRP qrp;
    qrp.circuit = circuit;
    qrp.midOffset = circuit.numberOfInputWires;
    qrp.outOffset = circuit.numberOfWires - circuit.numberOfOutputWires;
    qrp.t  = t;
    qrp.V = V;
    qrp.W = W;
    qrp.Y = Y;
    return qrp;
}

Vec<ZZ_p> eval(Circuit circuit, Vec<ZZ_p> input) {
    Vec<ZZ_p> allWireValues = Vec<ZZ_p>(input);
    allWireValues.SetLength(circuit.numberOfWires);
    for (int k = 0; k < circuit.numberOfMultiplicationGates; k++) {
        Gate kthGate = circuit.gates[k];
        ZZ_p leftInputSum;
        for (int i = 0; i < kthGate.leftInputs.size(); i++) {
            leftInputSum += allWireValues[kthGate.leftInputs[i]];
        }
        ZZ_p rightInputSum;
        for (int i = 0; i < kthGate.rightInputs.size(); i++) {
            rightInputSum += allWireValues[kthGate.rightInputs[i]];
        }
        allWireValues[k + circuit.numberOfInputWires] = leftInputSum * rightInputSum;
    }
    return allWireValues;
}