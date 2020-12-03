#include<qrp.h>
#include<gr.h>
#include<circuit.h>
#include<NTL/ZZ_pE.h>
#include<NTL/ZZX.h>
#include<NTL/tools.h>

using namespace NTL;
using namespace std;

Vec<ZZ_p> eval(const Circuit& circuit, const Vec<ZZ_p>& input) {
    Vec<ZZ_p> allWireValues = Vec<ZZ_p>(input);
    allWireValues.SetLength(circuit.numberOfWires);
    for (int k = 0; k < circuit.numberOfMultiplicationGates; k++) {
        // cout << "Evaluating wire " << k << endl;
        const vector<long> leftInputs = circuit.gates[k].leftInputs;
        const vector<long> rightInputs = circuit.gates[k].rightInputs;
        ZZ_p leftInputSum;
        for (long i = 0; i < leftInputs.size(); i++) {
            leftInputSum += allWireValues[leftInputs[i]];
        }
        ZZ_p rightInputSum;
        for (long i = 0; i < rightInputs.size(); i++) {
            rightInputSum += allWireValues[rightInputs[i]];
        }
        allWireValues[k + circuit.numberOfInputWires] = leftInputSum * rightInputSum;
    }
    return allWireValues;
}

void printCircuit(const Circuit& circuit) {
    cout << "numberOfWires..............: " << circuit.numberOfWires << "\n";
    cout << "numberOfInputWires.........: " << circuit.numberOfInputWires << "\n";
    cout << "numberOfMidWires...........: " << circuit.numberOfMidWires << "\n";
    cout << "numberOfOutputWires........: " << circuit.numberOfOutputWires << "\n";
    cout << "numberOfMultiplicationGates: " << circuit.numberOfMultiplicationGates << "\n";
    for (long i = 0; i < circuit.gates.size(); i++) {
        const vector<long> l = circuit.gates[i].leftInputs;
        const vector<long> r = circuit.gates[i].rightInputs;
        
        cout << "Gate " << i + circuit.numberOfInputWires << ":\n";
        cout << "  leftInputs:  ";
        for (long j = 0; j < l.size(); j++) {
            cout << " " << l[j];
        }
        cout << "\n";
        cout << "  rightInputs: ";
        for (long j = 0; j < r.size(); j++) {
            cout << " " << r[j];
        }
        cout << "\n";
        cout << "\n";
    }
}


ostream& operator<<(ostream& s, const Circuit& circuit) {
    s << circuit.numberOfWires << " ";
    s << circuit.numberOfInputWires << " ";
    s << circuit.numberOfOutputWires << "\n";
    for (long i = 0; i < circuit.gates.size(); i++) {
        const vector<long> l = circuit.gates[i].leftInputs;
        const vector<long> r = circuit.gates[i].rightInputs;
        s << l.size() << " " << r.size();
        
        for (long j = 0; j < l.size(); j++) {
            s << " " << l[j];
        }
        for (long j = 0; j < r.size(); j++) {
            s << " " << r[j];
        }
        s << "\n";
    }
    return s;
}

istream& operator>>(istream& s, Circuit& c) {
    if ((s >> c.numberOfWires) && 
        (s >> c.numberOfInputWires) && 
        (s >> c.numberOfOutputWires)) {
            c.numberOfMidWires = c.numberOfWires - c.numberOfInputWires - c.numberOfOutputWires;
            c.numberOfMultiplicationGates = c.numberOfMidWires + c.numberOfOutputWires;
    }

    long next;
    while (s >> next) {
        const long nlefts = next;
        if (!(s>> next)){
            return s; //todo error
        }
        const long nrights = next;
        vector<long> lefts, rights;
        for (long k = 0; k < nlefts; k++) {
            if (!(s >> next)) return s;//todo error
            lefts.push_back(next);
        }
        for (long k = 0; k < nrights; k++) {
            if (!(s >> next)) return s;//todo error
            rights.push_back(next);
        }
        c.gates.push_back(Gate{
            .leftInputs = lefts,
            .rightInputs = rights
        });
    }
    
    return s;
}
