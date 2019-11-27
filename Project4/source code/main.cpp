#include "isingmodel.hpp"
#include "omp.h"
//#include "mpi.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>


using namespace std;

int main() {


    double exp_MC_cyc = 10000.0; // 1 million
    bool to_file = false;

    cout << "Running OpenMP for parallelization!  "  << endl;

    cout << "The number of processors available = " << omp_get_num_procs()   << endl;
    cout << "The number of threads available    = " << omp_get_max_threads() <<  endl;
    cout << endl;

    // OpenMP used as problems with instaling MPI occured
#pragma omp parallel for
    for(int i = 200; i <= 230; i += 5) {

        double temp = i / 100.;

        cout << endl;
        cout << "Temperature =  " << temp  << " Thread number: " << omp_get_thread_num() << endl;

        IsingModel task_e_40 (40, "task_e_40");    task_e_40.initialize_random();
        IsingModel task_e_60 (60, "task_e_60");   task_e_60.initialize_random();
        IsingModel task_e_80 (80, "task_e_80");   task_e_80.initialize_random();
        IsingModel task_e_100(100, "task_e_100"); task_e_100.initialize_random();

        task_e_40.simulate(i, exp_MC_cyc, to_file);
        task_e_40.reset_expectations();

        task_e_60.simulate(i, exp_MC_cyc, to_file);
        task_e_60.reset_expectations();

        task_e_80.simulate(i, exp_MC_cyc, to_file);
        task_e_80.reset_expectations();

        task_e_100.simulate(i, exp_MC_cyc, to_file);
        task_e_100.reset_expectations();

    }



    return 0;
}
