#include<qrp.h>
#include<gr.h>
#include<NTL/ZZ_pE.h>
#include<NTL/ZZX.h>

#include<time.h>

using namespace NTL;
using namespace std;


ZZ_pX primitiveIrredPoly(long degree) {
    ZZ_pX P;
    P = 1;
    if (degree < 2) return P;
    if (degree > 32) return P;
    SetCoeff(P, degree);

    switch (degree)
    {
    case 2:
        // x^2 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 3:
        // x^3 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 4:
        // x^4 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 5:
        // x^5 + x^2 + 1
        SetCoeff(P, 2);
        break;
    case 6:
        // x^6 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 7:
        // x^7 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 8:
        // x^8 + x^4 + x^3 + x^2 + 1
        SetCoeff(P, 4);
        SetCoeff(P, 3);
        SetCoeff(P, 2);
        break;
    case 9:
        // x^9 + x^4 + 1
        SetCoeff(P, 4);
        break;
    case 10:
        // x^10 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 11:
        // x^11 + x^2 + 1
        SetCoeff(P, 2);
        break;
    case 12:
        // x^12 + x^6 + x^4 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 4);
        SetCoeff(P, 6);
        break;
    case 13:
        // x^13 + x^4 + x^3 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 3);
        SetCoeff(P, 4);
        break;
    case 14:
        // x^14 + x^8 + x^6 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 6);
        SetCoeff(P, 8);
        break;
    case 15:
        // x^15 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 16:
        // x^16 + x^9 + x^8 + x^7 + x^6 + x^4 + x^3 + x^2 + 1
        SetCoeff(P, 2);
        SetCoeff(P, 3);
        SetCoeff(P, 4);
        SetCoeff(P, 6);
        SetCoeff(P, 7);
        SetCoeff(P, 8);
        SetCoeff(P, 9);
        break;
    case 17:
        // x^17 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 18:
        // x^18 + x^5 + x^4 + x^3 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 3);
        SetCoeff(P, 4);
        SetCoeff(P, 5);
        break;
    case 19:
        // x^19 + x^5 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 5);
        break;
    case 20:
        // x^20 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 21:
        // x^21 + x^2 + 1
        SetCoeff(P, 2);
        break;
    case 22:
        // x^22 + x^1 + 1
        SetCoeff(P, 1);
        break;
    case 23:
        // x^23 + x^5 + 1
        SetCoeff(P, 5);
        break;
    case 24:
        // x^24 + x^7 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 7);
        break;
    case 25:
        // x^25 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 26:
        // x^26 + x^6 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 6);
        break;
    case 27:
        // x^27 + x^5 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 5);
        break;
    case 28:
        // x^28 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 29:
        // x^29 + x^2 + 1
        SetCoeff(P, 2);
        break;
    case 30:
        // x^30 + x^23 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 23);
        break;
    case 31:
        // x^31 + x^3 + 1
        SetCoeff(P, 3);
        break;
    case 32:
        // x^32 + x^22 + x^2 + x^1 + 1
        SetCoeff(P, 1);
        SetCoeff(P, 2);
        SetCoeff(P, 22);
        break;
    default:
        break;
    }
    return P;
}

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
    cout << "Subset of exceptional set computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";

    t = clock();
    const Vec<ZZ_pEX> termsOfT = getTargetPolynomialTerms(numberOfMultiplicationGates);
    cout << "terms of t computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n" ; 
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
    cout << "Target polynomial computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
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
    cout << "Deltas computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";//slow

    return deltas;
}

QRP getQRP(const Circuit& circuit, long k, long minDegree) {
    // ZZ mod = ZZ(1) << k;
    // ZZ_p::init(mod);
    // ZZ_pX polyMod = primitiveIrredPoly(minDegree); //todo check if minDegree is >= log_2(multGates)
    // ZZ_pE::init(polyMod);
    // todo fix automatic minimal modulus
    
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
    cout << "V, W and Y computed: " << ((double) t) / CLOCKS_PER_SEC << " seconds\n";
    
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
