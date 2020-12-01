#ifndef RINOCCHIO_H
#define RINOCCHIO_H

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>

#include<qrp.h>
#include<setup.h>
#include<joye_libert.h>

using namespace NTL;

//dummy encryption:
Vec<ZZ> E(ZZ_pE x);

// dummy decryption
ZZ_pE D(Vec<ZZ> y);

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

CRS getCRS(QRP prog, secretState ss);

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

Proof prove(QRP prog, CRS crs, Vec<ZZ_p> input);

bool verify(QRP qrp, secretState secret, CRS crs, Proof pi, Vec<ZZ_p> input, Vec<ZZ_p> output);



#endif