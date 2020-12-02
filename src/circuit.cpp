#include<qrp.h>
#include<gr.h>
#include<circuit.h>
#include<NTL/ZZ_pE.h>
#include<NTL/ZZX.h>

using namespace NTL;
using namespace std;

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