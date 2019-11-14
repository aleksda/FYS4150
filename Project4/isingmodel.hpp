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

        int num_of_spins;
        int E, M;

        Mat<int> spins;
        Col<int> index;

        vec boltzmann;	// Arma
        vec exp_value;	// Arma



        mt19937 rand_nums;
        uniform_real_distribution<double> uniform_zero_to_one;
        uniform_int_distribution<int> zero_to_L_distribution;

        IsingModel(int size, int data, string filename);

        //
        void initialize_random();
        void initialize_up();

        void get_lattice();
        void get_lattice_to_file(ofstream& file);

        double find_energy();
        double find_magnetization();

        int get_site_energy(int x, int y);

        double metropolis();

        void equilibrate(double temp, int num_of_mc_cycles);
        void simulate(double temp, int num_of_mc_cycles, bool to_file);

        void reset_expectations();

        void simulate_to_file(double temp);
};

#endif // ISINGMODEL_H
