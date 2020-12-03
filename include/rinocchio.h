#ifndef RINOCCHIO_H
#define RINOCCHIO_H

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>

#include<qrp.h>
#include<setup.h>
#include<joye_libert.h>

using namespace NTL;

struct CRS {
    Vec<JLEncoding> powersOfS;
    Vec<JLEncoding> powersOfSMultAlpha;
    Vec<JLEncoding> rvVofS;
    Vec<JLEncoding> rwWofS;
    Vec<JLEncoding> ryYofS;
    Vec<JLEncoding> alpharvVofS;
    Vec<JLEncoding> alpharwWofS;
    Vec<JLEncoding> alpharyYofS;
    Vec<JLEncoding> betaSums; 
    JLEncodingKey publicKey; //todo separate public/secret
};

CRS getCRS(const QRP& prog, const secretState& ss);

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

Proof prove(const QRP& prog, const CRS& crs, const Vec<ZZ_p>& allWireValues);

bool verify(const QRP& qrp, const secretState& secret, const CRS& crs, const Proof& pi, const Vec<ZZ_p>& input, const Vec<ZZ_p>& output);



#endif