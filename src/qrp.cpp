#include<qrp.h>
#include<gr.h>
#include<NTL/ZZ_pE.h>
#include<NTL/ZZX.h>

using namespace NTL;
using namespace std;

QRP getQRP(Circuit circuit) {
    //Pick distinct elements of exceptional set for each gate:
    Vec<ZZ_pE> multWires;


    // Compute t(x) = (x - g_1) * (x - g_2)
    Vec<ZZ_pEX> V, W, Y;
    V.SetLength(circuit.numberOfWires);
    W.SetLength(circuit.numberOfWires);
    Y.SetLength(circuit.numberOfWires);
    ZZ_pEX t = conv<ZZ_pEX>(1);
    ZZ_pE galloisOne;
    set(galloisOne);
    ZZ_pEX x;
    SetX(x);
    
    Vec<Vec<ZZ_pE>> v_values, w_values, y_values;
    v_values.SetLength(circuit.numberOfWires);
    w_values.SetLength(circuit.numberOfWires);
    y_values.SetLength(circuit.numberOfWires);
    
    for (int k = 0; k < circuit.numberOfMultiplicationGates; k++) {
        multWires.append(indexedElementInExceptionalSet(k));

        //target polynomial:
        t *= x - indexedElementInExceptionalSet(k);

        Gate kthMultGate =  circuit.gates[k];

        int nextLeftInput = 0;
        int nextRightInput = 0;
        vector<long> leftInputs = kthMultGate.leftInputs;
        vector<long> rightInputs = kthMultGate.rightInputs;
        for (int j = 0; j < circuit.numberOfWires; j++) {
            if (j == k + circuit.numberOfInputWires) {
                y_values[j].append(galloisOne);
            } else {
                y_values[j].append(ZZ_pE::zero());
            }
            
            if (nextLeftInput < leftInputs.size() && j == leftInputs[nextLeftInput]) {
                v_values[j].append(galloisOne);
                nextLeftInput++;
            } else {
                v_values[j].append(ZZ_pE::zero());
            }
            
            if (nextRightInput < rightInputs.size() && j == rightInputs[nextRightInput]) {
                w_values[j].append(galloisOne);
                nextRightInput++;
            } else {
                w_values[j].append(ZZ_pE::zero());
            }
        }
    }

    for (int k = 0; k < circuit.numberOfWires; k++) {
        interpolate(V[k], multWires, v_values[k]);
        interpolate(W[k], multWires, w_values[k]);
        interpolate(Y[k], multWires, y_values[k]);
        
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