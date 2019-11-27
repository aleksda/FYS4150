#ifndef ISINGMODEL_H
#define ISINGMODEL_H

#include <armadillo>
#include <random>
#include <string>

using namespace std; using namespace arma;


class IsingModel {

	private:

	    int get_site_energy(int x, int y);

    public:

        int size;
        //int data;
        string filename;

        int num_of_spins;
        int E, M;

        Mat<int> spins;
        Col<int> index;

        vec boltzmann;	// Arma
        vec exp_mag_etc_values;	// Arma

        // Needed because rand() did not work (as in lecture notes)
        mt19937 ran_vals;
        uniform_real_distribution<double> uniform_zero_to_one;
        uniform_int_distribution<int> rand_num;

        IsingModel(int size, string filename);

        void initialize_random();
        void initialize_up();

        void get_lattice();
        void get_lattice_to_file(ofstream& file);

        double find_energy();
        double find_magnetization();

        void metropolis();

        void equilibrate(double temp, int num_of_mc_cycles);
        void simulate(double temp, int num_of_mc_cycles, bool to_file);

        void reset_expectations();

        //void simulate_to_file(double temp);
        void simulate_to_file(int num_of_mc_cycles, ofstream &ofile_except);
        void time_to_file(bool start);
};

#endif // ISINGMODEL_H
