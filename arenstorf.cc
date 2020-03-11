//
// Created by andie on 09.03.2020.
//

#include "arenstorf.h"
#include "DormandPrinceIntegrator.h"
#include <iostream>
#include <time.h>

void
Model_Arenstorf::run() {
    long double t_start = 0;
    long double period = 11.124340337266085134999734047;
    long double t_finish = 1;
    long double delta = 0.01;
    long double time;
    int size = (int)((t_finish - t_start) / delta) +1;
    result_list = new long double* [size];
    for (int i = 0; i < size; i++) {
        result_list[i] = new long double[n + 1];
    }
    Integrator_Dormand_Prince integrator(1e-16, 1e-16, 1e-16, 1e-16, n, 0.01, size);
    time = clock();
    long double* result = integrator.integration(t_start, t_finish, this);
    time = (clock() - time) / CLOCKS_PER_SEC;
    std::cout << "time elapsed " << time << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << result[i] << std::endl;
    }
    plotter(11, 1, size);


    integrator.~Integrator_Dormand_Prince();
}

void
Model_Arenstorf::function(long double* ind, long double* dep,
        long double* step_result) {
 // 1 in right parts ~R0;
    Ic= k_B*CalculateCoefficientIb(dep[1]);
    step_result[0] =(dep[2]-Ic-dep[3])/C1;
    //std::cout<<"st_res 0 "<<step_result[7]<<std::endl;
    step_result[1]=((Ve-dep[1])/Re -dep[2]-CalculateCoefficientIb(dep[1]))/C2;
    step_result[2]=(Vc-dep[0]-Rl*dep[2]+dep[1])/L;
    step_result[3]=dep[4];
    step_result[4]=((dep[2]+Ic)/C1+dep[5]/C0-R0*dep[4]-(1/C+1/C0+1/C1)*dep[3])/L0;
    step_result[5]=dep[6];
    step_result[6]=((dep[3]+dep[7])/C0-R0*dep[6]-(1/C+2/C0)*dep[5])/L0;
    step_result[7]=dep[8];
    step_result[8]=((dep[5]+dep[9])/C0-(1/C+2/C0)*dep[7]-R0*dep[8])/L0;
    step_result[9]=dep[10];
    step_result[10]=(dep[7]/C0-(1/C+2/C0)*dep[9]-R0*dep[10])/L0;

}

long double
Model_Arenstorf::CalculateCoefficientIb( long double dep) {
    long double result=0.0;
    if (dep<=Vthr){
        result=0;

        return result;
    }
    else{
        result=alpha*(dep-Vthr);

        return result;
    }
}

void
Model_Arenstorf::spectralPower() {

}
