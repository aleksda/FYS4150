#include "isingmodel.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>


using namespace std;

int main() {


    double exp_MC_cyc = 32.0;
    bool to_file = true;
/*
    int L_b = 2;
    double T_b = 1.0;

    IsingModel test_b(L_b, 10000, "test_b");

    test_b.initialize_up();
    test_b.get_lattice();

    test_b.equilibrate(T_b, exp_MC_cyc);
    test_b.simulate(T_b, pow(10, exp_MC_cyc), to_file);
    test_b.reset_expectations();
    cout << endl;
    test_b.get_lattice();

    for (double t = 0.5; t < 4.0; t += 0.05) {
      test_b.equilibrate(t, exp_MC_cyc);
      test_b.simulate(t, pow(10, exp_MC_cyc), to_file);
      test_b.reset_expectations();
    }
*/

    int L_c = 20;

    double T_c1 = 1.0;


    IsingModel task_c_unordered(L_c, pow(10, exp_MC_cyc - 2), "task_c_unordered"); // -2 ?
/*

    task_c_unordered.initialize_random();
    task_c_unordered.get_lattice();

    // Test with T = 1.0
    task_c_unordered.initialize_random();
    task_c_unordered.equilibrate(T_c1, exp_MC_cyc);
    task_c_unordered.simulate(T_c1, pow(10, exp_MC_cyc), to_file);
    task_c_unordered.reset_expectations();
*/

    double T_c2 = 2.4;

    task_c_unordered.initialize_random();
    task_c_unordered.get_lattice();

    // Test with T = 2.4
    task_c_unordered.initialize_random();
    task_c_unordered.equilibrate(T_c2, exp_MC_cyc);
    task_c_unordered.simulate(T_c1, pow(10, exp_MC_cyc), to_file);
    task_c_unordered.reset_expectations();

/*
    IsingModel task_c_ordered(L_c, "task_c_ordered", pow(10, exp_MC_cyc - 1)); // -1 ?

    // Test with T = 1.0
    task_c_ordered.initialize_up();
    task_c_ordered.equilibrate(T_c1, exp_MC_cyc);
    task_c_ordered.simulate(T_c1, pow(10, exp_MC_cyc), to_file);
    task_c_ordered.reset_expectations();

    // Test with T = 2.4
    task_c_unordered.initialize_random();
    task_c_unordered.equilibrate(T_c2, exp_MC_cyc);
    task_c_unordered.simulate(T_c2, pow(10, exp_MC_cyc), to_file);
    task_c_unordered.reset_expectations();

*/

    return 0;
}
