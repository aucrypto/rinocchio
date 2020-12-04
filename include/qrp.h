#ifndef QRP_H
#define QRP_H

#include <NTL/ZZ_pEX.h>
#include <NTL/vector.h>
#include <circuit.h>

using namespace NTL;




struct QRP {
    Circuit circuit;
    long midOffset, outOffset;
    ZZ_pEX t; //target polynomial
    Vec<ZZ_pEX> V;
    Vec<ZZ_pEX> W;
    Vec<ZZ_pEX> Y;
};

QRP getQRP(const Circuit& circuit);

ostream& operator<<(ostream& s, const QRP& qrp);
istream& operator>>(istream& s, QRP& qrp);

#endif