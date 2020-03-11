#include <iostream>
#include "Object.h"
#include "DormandPrinceIntegrator.h"
#include "arenstorf.h"
int
main() {




        long double *start_conditions = new long double[11];
        for ( int i =0 ; i<11; ++i ){
            start_conditions[i]=0;
        }

        Model_Arenstorf *b = new Model_Arenstorf(11, start_conditions);
        b->run();
    }



