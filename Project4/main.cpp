#include "isingmodel.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>


using namespace std;

int main() {

    int L = 2;

    double exp_MC_cyc = 32.0;

    bool to_file = true;

    IsingModel test_b(L, 10000, "test_b");
    test_b.initialize_up();
    test_b.get_lattice();
/*
    for (double t = 0.5; t < 4.0; t += 0.05) {
      test_b.equilibrate(t, exp_MC_cyc);
      test_b.simulate(t, pow(10, exp_MC_cyc), to_file);
      test_b.reset_expectations();
    }
*/
    //ofstream test("test.txt");
    //if (!test.is_open()) cout << "Could not open test file" << endl;
    ofstream test;
    test.open("test.txt");
    test << " test passed!" << "\n";
    test.close();

    return 0;
}
