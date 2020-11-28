#include<qrp.h>
#include<gr.h>
#include<NTL/ZZ_pE.h>
#include<NTL/ZZX.h>

using namespace NTL;
using namespace std;

QRP getQRP(Circuit circuit) {
    //Pick distinct elements of exceptional set for each gate:
    ZZ_pE g_1, g_2;
    g_1 = indexedElementInExceptionalSet(1);
    g_2 = indexedElementInExceptionalSet(2);
    ZZ_pE diff = g_1 - g_2;

    // Compute t(x) = (x - g_1) * (x - g_2)
    Vec<ZZ_pEX> V, W, Y;
    V.SetLength(6);
    W.SetLength(6);
    Y.SetLength(6);
    ZZ_pEX t;
    
    {
        ZZ_pEX x;
        SetX(x);
        t = (x - g_1)*(x-g_2);

        ZZ_pE galloisOne;
        set(galloisOne);

        Vec<ZZ_pE> a;
        a.append(g_1);
        a.append(g_2);
        Vec<ZZ_pE> b_v_1;
        b_v_1.append(ZZ_pE::zero());
        b_v_1.append(galloisOne);
        Vec<ZZ_pE> b_v_2;
        b_v_2.append(ZZ_pE::zero());
        b_v_2.append(galloisOne);
        Vec<ZZ_pE> b_v_3;
        b_v_3.append(galloisOne);
        b_v_3.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_v_4;
        b_v_4.append(ZZ_pE::zero());
        b_v_4.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_v_5;
        b_v_5.append(ZZ_pE::zero());
        b_v_5.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_v_6;
        b_v_6.append(ZZ_pE::zero());
        b_v_6.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_w_1;
        b_w_1.append(ZZ_pE::zero());
        b_w_1.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_w_2;
        b_w_2.append(ZZ_pE::zero());
        b_w_2.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_w_3;
        b_w_3.append(ZZ_pE::zero());
        b_w_3.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_w_4;
        b_w_4.append(galloisOne);
        b_w_4.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_w_5;
        b_w_5.append(ZZ_pE::zero());
        b_w_5.append(galloisOne);
        Vec<ZZ_pE> b_w_6;
        b_w_6.append(ZZ_pE::zero());
        b_w_6.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_1;
        b_y_1.append(ZZ_pE::zero());
        b_y_1.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_2;
        b_y_2.append(ZZ_pE::zero());
        b_y_2.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_3;
        b_y_3.append(ZZ_pE::zero());
        b_y_3.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_4;
        b_y_4.append(ZZ_pE::zero());
        b_y_4.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_5;
        b_y_5.append(galloisOne);
        b_y_5.append(ZZ_pE::zero());
        Vec<ZZ_pE> b_y_6;
        b_y_6.append(ZZ_pE::zero());
        b_y_6.append(galloisOne);

        interpolate(V(1), a, b_v_1);
        interpolate(V(2), a, b_v_2);
        interpolate(V(3), a, b_v_3);
        interpolate(V(4), a, b_v_4);
        interpolate(V(5), a, b_v_5);
        interpolate(V(6), a, b_v_6);

        interpolate(W(1), a, b_w_1);
        interpolate(W(2), a, b_w_2);
        interpolate(W(3), a, b_w_3);
        interpolate(W(4), a, b_w_4);
        interpolate(W(5), a, b_w_5);
        interpolate(W(6), a, b_w_6);

        interpolate(Y(1), a, b_y_1);
        interpolate(Y(2), a, b_y_2);
        interpolate(Y(3), a, b_y_3);
        interpolate(Y(4), a, b_y_4);
        interpolate(Y(5), a, b_y_5);
        interpolate(Y(6), a, b_y_6);
    }
    
    QRP qrp;
    qrp.circuit = circuit;
    qrp.midOffset = circuit.numberOfInputWires;
    qrp.outOffset = circuit.numberOfWires - circuit.numberOfOutputWires;
    qrp.t  = t;
    qrp.V = V;
    qrp.W = W;
    qrp.Y = Y;
    return qrp;
}