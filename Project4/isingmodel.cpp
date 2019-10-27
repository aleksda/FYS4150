#include "isingmodel.hpp"

#include <armadillo>
#include <random>
#include <string>

using namespace std; using namespace arma;

IsingModel::IsingModel(int size, int data, string filename) {
    // L = size = latice_dim

    this->size     = size;
    this->data     = data;
    this->filename = filename + ".csv";

    int num_of_spins   = size * size;

    spins = Mat<int>(size, size);

    // From std library
    mt19937 rand_nums(clock());
    uniform_real_distribution<double> uniform_zero_to_one(0.0, 1.0);
    //stduniform_int_distribution<int> uniform_zero_to_one(0.0, 1.0);


}

void IsingModel::initialize_random() {

    int pm_one = 0;

    for(unsigned int i = 0; i < this->size; i++) {
        for(unsigned int j = 0; j < this->size; j++) {
            if(uniform_zero_to_one(rand_nums) >= 0.5) pm_one = 1;
            else pm_one = -1;
            spins(i, j) = pm_one;
        }
    }
}

// Not part of project
void IsingModel::get_lattice() {
    //for(unsigned int i = 0; i < this->size; i++) {
        //7for(unsigned int j = 0; j < this->size; j++) {

            spins.print(cout);
        //}
    //}
}
