#ifndef RINOCCHIO_H
#define RINOCCHIO_H

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>

#include<qrp.h>
#include<setup.h>
#include<joye_libert.h>

using namespace NTL;

struct Proof {
    JLEncoding rvVmidOfS;
    JLEncoding rwWmidOfS;
    JLEncoding ryYmidOfS;
    JLEncoding alphaVrvVmidOfS;
    JLEncoding alphaWrwWmidOfS;
    JLEncoding alphaYryYmidOfS;
    JLEncoding betaSum;
    JLEncoding hOfS;
    JLEncoding alphaHOfS;
};

ZZ_pEX proverComputeH(const QRP& prog, const Vec<ZZ_p>& allWireValues);

Proof prove(const ZZ_pEX& H, const CRS& crs, const Vec<ZZ_p>& allWireValues, long midOffset, long outOffset);

Proof prove(const QRP& prog, const CRS& crs, const Vec<ZZ_p>& allWireValues);

bool verify(const SecretState& secret, const CRS& crs, const Proof& pi, const Vec<ZZ_p>& input, const Vec<ZZ_p>& output);



#endif