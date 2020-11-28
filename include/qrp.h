#ifndef QRP_H
#define QRP_H

#include <NTL/ZZ_pEX.h>

struct Circuit {
    long numberOfWires;
    long numberOfInputWires;
    long numberOfMidWires;
    long numberOfOutputWires;
    long numberOfMultiplicationGates;
    
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