#ifndef SETUP_H
#define SETUP_H

#include <NTL/ZZ_pE.h>
#include <joye_libert.h>
#include <qrp.h>

using namespace NTL;

struct SecretState {
    ZZ_pE s, r_v, r_w, r_y, alpha, alpha_v, alpha_w, alpha_y, beta;
    Vec<ZZ_pE> inVofS, inWofS, inYofS, outVofS, outWofS, outYofS;
    ZZ_pE tOfS;
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
