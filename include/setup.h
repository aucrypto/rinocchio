#ifndef SETUP_H
#define SETUP_H

#include <NTL/ZZ_pE.h>
#include <joye_libert.h>
#include <qrp.h>

using namespace NTL;

struct SecretState {
    ZZ_pE s;
    ZZ_pE r_v;
    ZZ_pE r_w;
    ZZ_pE r_y;
    ZZ_pE alpha;
    ZZ_pE alpha_v;
    ZZ_pE alpha_w;
    ZZ_pE alpha_y;
    ZZ_pE beta;
    ZZ_pE tOfS;
    Vec<ZZ_pE> inVofS;
    Vec<ZZ_pE> inWofS;
    Vec<ZZ_pE> outVofS;
    Vec<ZZ_pE> outWofS;
    Vec<ZZ_pE> outYofS;
    JLEncodingKey secretKey;
};

ostream& operator<<(ostream& s, const SecretState& secretState);
istream& operator>>(istream& s, SecretState& secretState);

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
    JLEncodingKey publicKey;
};


ostream& operator<<(ostream& s, const CRS& crs);
istream& operator>>(istream& s, CRS& crs);

SecretState setup(const QRP& qrp, long l, long k);

CRS getCRS(const QRP& qrp, const SecretState& secret);


#endif
