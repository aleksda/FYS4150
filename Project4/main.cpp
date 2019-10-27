#include "isingmodel.hpp"

#include <iostream>

using namespace std;

int main() {

    int L = 2;

    IsingModel test_b(L, NULL, "test_b");
    test_b.initialize_random();
    test_b.get_lattice();

    return 0;
}
