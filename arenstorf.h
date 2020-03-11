//
// Created by andie on 09.03.2020.
//

#ifndef KR_ARENSTORF_H
#define KR_ARENSTORF_H

#include "Object.h"

class Model_Arenstorf : public Object {
public:

    Model_Arenstorf(int n_, long double*
    start_conditions)
            : Object(11, start_conditions) {
        Vthr = 0.75;
        C1 = C2 = C0 = 70e-5;
        C = 140e-5;
        L = 1;
        alpha = 1 / 70;
        Vc = Ve = 20;
        Rl = 95;
        Re = 2.4*1000;
        L0 = 5e-6;
        k_B = 20;
        R0 = 0.00001;
        Ic = 1;

    };

    virtual void function(long double* ind, long double* dep,
            long double* step_result) override;

    virtual long double* GetA() override { return a; }

    virtual long double GetA_I(int i) override { return a[i]; }

    void spectralPower();
    virtual void run();
    ~Model_Arenstorf() = default;
protected:
    long double L;
    long double k_B;
    long double C1, C2, C0, C;
    long double Vthr;
    long double alpha, Ic, Ib, Rl, Re, Vc, Ve;
    long double L0, R0;

    long double
    CalculateCoefficientIb(long double dep);
};

#endif //KR_ARENSTORF_H
