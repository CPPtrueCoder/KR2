//
// Created by andie on 08.03.2020.
//

#ifndef KR_OBJECT_H
#define KR_OBJECT_H

#include <iostream>
#include <cmath>

class Object {
protected:
    long double *a;
    int n;

public:
    Object() = default;
    Object(int n_, long double *start_conditions);
    ~Object();
    virtual void run(long double t_start, long double t_finish);
    virtual void function(long double *ind, long double *dep,
            long double *step_result);
    virtual long double *GetA() { return a; }
    virtual void SetA(long double *a_);
    virtual long double GetA_I(int i) { return a[i]; }
    long double **result_list = nullptr;
    void plotter(int k, int l, int size);
};

class NewObject : public Object {

    void run(long double t_start, long double t_finish) override;
};

#endif //KR_OBJECT_H
