#ifndef RINOCCHIO_H
#define RINOCCHIO_H

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>

#include<qrp.h>
#include<setup.h>

using namespace NTL;

//dummy encryption:
Vec<ZZ> E(ZZ_pE x);

// dummy decryption
ZZ_pE D(Vec<ZZ> y);

struct CRS {
    Vec<Vec<ZZ>> powersOfS;
    Vec<Vec<ZZ>> powersOfSMultAlpha;
    Vec<Vec<ZZ>> rvVofS;
    Vec<Vec<ZZ>> rwWofS;
    Vec<Vec<ZZ>> ryYofS;
    Vec<Vec<ZZ>> alpharvVofS;
    Vec<Vec<ZZ>> alpharwWofS;
    Vec<Vec<ZZ>> alpharyYofS;
    Vec<Vec<ZZ>> betaSums; 
    bool publicKey;
};

CRS getCRS(QRP prog, secretState ss);

struct Proof {
    Vec<ZZ> rvVmidOfS;
    Vec<ZZ> rwWmidOfS;
    Vec<ZZ> ryYmidOfS;
    Vec<ZZ> alphaVrvVmidOfS;
    Vec<ZZ> alphaWrwWmidOfS;
    Vec<ZZ> alphaYryYmidOfS;
    Vec<ZZ> betaSum;
    Vec<ZZ> hOfS;
    Vec<ZZ> alphaHOfS;
};

Proof prove(QRP prog, CRS crs, Vec<ZZ_p> input);

bool verify(QRP qrp, secretState secret, CRS crs, Proof pi, Vec<ZZ_p> input, Vec<ZZ_p> output);



#endif