#ifndef SETUP_H
#define SETUP_H

#include <NTL/ZZ_pE.h>
#include <joye_libert.h>

using namespace NTL;

struct secretState {
    ZZ_pE s, r_v, r_w, r_y, alpha, alpha_v, alpha_w, alpha_y, beta;
    JLEncodingKey secretKey;
};

secretState setup(long l, long k);

#endif
