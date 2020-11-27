#ifndef QRP_H
#define QRP_H

#include <NTL/ZZ_pEX.h>

struct QRP {
    long numberOfWires;
    long numberOfInputWires;
    long numberOfOutputWires;
    NTL::ZZ_pEX t; //target polynomial
    NTL::Vec<NTL::ZZ_pEX> V;
    NTL::Vec<NTL::ZZ_pEX> W;
    NTL::Vec<NTL::ZZ_pEX> Y;
};

QRP getQRP(); //todo should take a circuit somehow

#endif