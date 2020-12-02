#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <NTL/ZZ_pEX.h>
#include <NTL/vector.h>
#include <vector>

using namespace NTL;

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

Vec<ZZ_p> eval(Circuit circuit, Vec<ZZ_p> input);

#endif
