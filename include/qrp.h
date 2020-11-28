#ifndef QRP_H
#define QRP_H

#include <NTL/ZZ_pEX.h>
#include <vector>

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

struct QRP {
    Circuit circuit;
    long midOffset, outOffset;
    NTL::ZZ_pEX t; //target polynomial
    NTL::Vec<NTL::ZZ_pEX> V;
    NTL::Vec<NTL::ZZ_pEX> W;
    NTL::Vec<NTL::ZZ_pEX> Y;
};

QRP getQRP(Circuit circuit);

#endif