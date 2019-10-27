#ifndef ISINGMODEL_H
#define ISINGMODEL_H

#include <armadillo>
#include <random>
#include <string>

using namespace std; using namespace arma;


class IsingModel {

    public:

        int size;
        int data;
        string filename;

        Mat<int> spins;

        mt19937 rand_nums;
        uniform_real_distribution<double> uniform_zero_to_one;

        IsingModel();
        IsingModel(int size, int data, string filename);

        //
        void initialize_random();
        void get_lattice();
};

#endif // ISINGMODEL_H
