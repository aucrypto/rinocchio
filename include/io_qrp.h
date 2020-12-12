#ifndef IOQRP_H
#define IOQRP_H

#include <NTL/ZZ_pEX.h>
#include <NTL/vector.h>
#include <circuit.h>

using namespace NTL;




struct IOQRP {
    Circuit circuit;
    long midOffset, outOffset;
    ZZ_pEX t; //target polynomial
};

IOQRP writeIOQRP(string path, const Circuit& circuit, long k, long minimumDegree);
void writeYPolys(string path, const IOQRP& qrp);
void writeVandWPolys(string vPath, string wPath, string yPath, const IOQRP& qrp);

ostream& operator<<(ostream& s, const IOQRP& qrp);
istream& operator>>(istream& s, IOQRP& qrp);

#endif