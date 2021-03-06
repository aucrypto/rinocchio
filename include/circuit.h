#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <NTL/ZZ_pEX.h>
#include <NTL/vector.h>
#include <vector>

using namespace NTL;
using namespace std;

struct Gate {
    std::vector<long> leftInputs, rightInputs;
};

struct Circuit {
    long numberOfWires;
    long numberOfInputWires;
    long numberOfMidWires;
    long numberOfOutputWires;
    long numberOfMultiplicationGates;
    std::vector<Gate> gates;
};

void printCircuit(const Circuit& circuit);

ostream& operator<<(ostream& s, const Circuit& circuit);
istream& operator>>(istream& s, Circuit& circuit);

Vec<ZZ_p> eval(const Circuit& circuit, const Vec<ZZ_p>& input);

#endif
