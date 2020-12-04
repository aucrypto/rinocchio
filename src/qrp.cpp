#include<qrp.h>
#include<gr.h>
#include<NTL/ZZ_pE.h>
#include<NTL/ZZX.h>

#include<time.h>

using namespace NTL;
using namespace std;

Vec<ZZ_pEX> getInterpolationDeltas(long numberOfMultiplicationGates, ZZ_pEX& targetPolynomial) {
    Vec<ZZ_pEX> deltas;
    deltas.SetLength(numberOfMultiplicationGates);

    ZZ_pEX x;
    SetX(x);
    set(targetPolynomial);

    cout << "Computing subset of exceptional set...\n";
    clock_t t;
    const Vec<ZZ_pE> exceptionalSubSet = getExceptionalSubset(numberOfMultiplicationGates);
    t = clock() - t;
    cout << "Subset of exceptional set computed: " << t / CLOCKS_PER_SEC << " seconds\n";

    t = clock();
    const Vec<ZZ_pEX> termsOfT = getTargetPolynomialTerms(numberOfMultiplicationGates);
    cout << "terms of t computed: " << t / CLOCKS_PER_SEC << " seconds\n" ; 
    //We get here very fast

    // Vec<ZZ_pE> exceptionalSubSet;
    // exceptionalSubSet.SetLength(numberOfMultiplicationGates);
    t = clock();
    for (long i = 0; i < numberOfMultiplicationGates; i++)  {
        // cout << "Compute exceptional elements and target: " << i << endl;
        // const ZZ_pE exceptionalElement = exceptionalSubSet[i];
        // targetPolynomial *= x - exceptionalElement;
        targetPolynomial *= termsOfT[i];//todo I think this was slower
    }
    t = clock() - t;
    cout << "Target polynomial computed: " << t / CLOCKS_PER_SEC << " seconds\n";
    //slow, would using a c style array be any faster?

    t = clock();
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
    t = clock() - t;
    cout << "Deltas computed: " << t / CLOCKS_PER_SEC << " seconds\n";//slow

    return deltas;
}

QRP getQRP(const Circuit& circuit) {
    //Pick distinct elements of exceptional set for each gate:
    ZZ_pEX targetPolynomial;
    Vec<ZZ_pEX> deltas = getInterpolationDeltas(circuit.numberOfMultiplicationGates, targetPolynomial);

    clock_t t;
    Vec<ZZ_pEX> V, W, Y;
    V.SetLength(circuit.numberOfWires);
    W.SetLength(circuit.numberOfWires);
    Y.SetLength(circuit.numberOfWires);

    ZZ_pEX x;
    SetX(x);
    t = clock();
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
    t = clock() - t;
    cout << "V, W and Y computed: " << t / CLOCKS_PER_SEC << " seconds\n";
    
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
