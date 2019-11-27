#include "isingmodel.hpp"

#include <armadillo>
#include <iomanip>
#include <random>
#include <string>
#include <chrono>

using namespace std; using namespace arma;
using namespace std::chrono;

stringstream stream;

/*
ofstream ofile_spin;   ofstream oflie_expect;
*/

IsingModel::IsingModel(int size, string filename) {

    this->size     = size;
    this->filename = filename; // + ".csv"

    int num_of_spins = size * size;
    int E, M = 0;

    // vectors & matrices
    spins = Mat<int>(size, size); // Mat<> before spins not working?
    index = Col<int>(size + 2);

    boltzmann  = vec(17, fill::zeros);
    exp_mag_etc_values  = vec(5, fill::zeros);

    // Random number ran_valss
    // From std library
    mt19937 ran_vals(clock()); // CHANGE to high precission clock ?
    uniform_real_distribution<double> uniform_zero_to_one(0.0, 1.0);
    uniform_int_distribution<int> rand_num(0, size-1); // rand() not correct

    index(0) = size - 1;
    for(int i = 1; i < size + 2; i++)
        index(i) = i - 1;

    index(size + 1) = 0;
}

void IsingModel::initialize_random() {

    //spins.randu()(); // Not working?
    int pm_one = 0;

    for(unsigned int i = 0; i < this->size; i++) {
        for(unsigned int j = 0; j < this->size; j++) {
            if(uniform_zero_to_one(ran_vals) >= 0.5) pm_one = 1;
            else pm_one = -1;
            spins(i, j) = pm_one;
        }
    }
}

void IsingModel::initialize_up() {
    spins.ones();
}

void IsingModel::get_lattice() {
    //for(unsigned int i = 0; i < this->size; i++) {
        //for(unsigned int j = 0; j < this->size; j++) {

            spins.print(cout);
        //}
    //}
}

void IsingModel::get_lattice_to_file(ofstream& file) {
    spins.save(file, arma_ascii);
}

// Calculates the total energy of given system
double IsingModel::find_energy() {
    int energy = 0;
    for (unsigned int i = 0; i < this->size; i++) {
        for (unsigned int j= 0; j < this->size; j++) {
            energy -= spins(i, j) * (spins(index(i+1), index(j)) \
                + spins(index(i), index(j+1)));
        }
    }

    return energy;
}

// Calculates the total magnetization of given system
double IsingModel::find_magnetization() {
    int magnetization = 0;
    for (unsigned int i = 0; i < this->size; i++) {
        for (unsigned int j = 0; j < this->size; j++)
            magnetization += spins(i, j);
    }

    return magnetization;
}

int IsingModel::get_site_energy(int x, int y) {

    int up 		 = spins(index(y), index(x+1));
    int down 	 = spins(index(y), index(x-1));
    int right  	 = spins(index(y+1), index(x));
    int left 	 = spins(index(y-1), index(x));
    int current  = spins(index(y), index(x));

    return 2 * current * (up + down + right + left);
}

// Metropolis Algorithm
void IsingModel::metropolis() {
    /*

    // This does not work? rand() not correct ?
    int delta_E;
    int rand_x; int rand_y;

    for(int i = 0; i < size * size; i++) {
        rand_x = rand()%(size+1); rand_y = rand()%(size+1);

        delta_E = get_site_energy(rand_y, rand_x);

        if(delta_E < 0 || (double) rand() <= boltzmann(delta_E + 8)) {

            E += delta_E;
            M -= 2 * spins(index(rand_y), index(rand_x));

            spins(index(rand_y), index(rand_x)) *= -1;
        }

    }
    */

    int delta_E;

    for (int i = 0; i < size * size; i++) {

        int rand_x = rand_num(ran_vals)%size+1;
        int rand_y = rand_num(ran_vals)%size+1;

        int get_site_energy = spins(index(rand_y), index(rand_x + 1)) + spins(index(rand_y), index(rand_x - 1)) + spins(index(rand_y+1), index(rand_x)) + spins(index(rand_y-1), index(rand_x));

        delta_E = 2 * get_site_energy * spins(index(rand_y), index(rand_x));

        if (delta_E < 0 || uniform_zero_to_one(ran_vals) <= boltzmann(delta_E + 8)) {

            E += delta_E;
            M -= 2 * spins(index(rand_y), index(rand_x));

            spins(index(rand_y) , index(rand_x)) *= -1;
        }
    }
}

void IsingModel::equilibrate(double temp, int num_of_mc_cycles) { // Better name?

    double beta = 1.0 / temp;

    boltzmann(0)  = 	exp(8.0 * beta);	// -8 +8
    boltzmann(4)  = 	exp(4.0 * beta);	// -4 +8
    boltzmann(8)  = 	1;					//  0 +8
    boltzmann(12) = 	exp(-4.0 * beta);	//  4 +8
    boltzmann(16) = 	exp(-8.0 * beta);	//  8 +8

    /*
    for(int i = 0; i <= 16; i += 4) {
        boltznann(i) = exp(8.0 * beta);
    }
    */

    for (int i = 0; i < num_of_mc_cycles / 10; i++) {
        metropolis();
    }
}

void IsingModel::simulate(double temp, int num_of_mc_cycles, bool to_file) {

    E = find_energy(); M = find_magnetization();

    // Data used to be a variable given in the constructor. Now changed
    //int data = 1000; // change?

    //int denominator = num_of_mc_cycles / data; //this->data;
    double beta = 1.0 / temp;
/*
    if(to_file)
        simulate_to_file(temp);
*/
    //if(to_file) {
        stream << fixed << setprecision(3) << temp;

        //ofile_spin.open("data/" + this->filename + "_" + to_string(stream) + "_spin.txt");
        ofstream ofile_spin("../data/" + this->filename + "_" + stream.str() + "_spin.txt");
        if (!ofile_spin.is_open()) cout << "Directory not found. Make sure you are in correct directory! " << endl;

        //oflie_expect.open("data/" + this->filename + "_" + to_string(stream) + "_expect.txt");
        ofstream oflie_expect("../data/" + this->filename + "_" + stream.str() + "_expectations.txt");
        if (!oflie_expect.is_open()) cout << "Directory not found. Make sure you are in correct directory! " << endl;

        get_lattice_to_file(ofile_spin);
    //}

    // Trick for reducing the FLOPs a little during computation
    int size_squared 		 = size * size;	// Alternativly use num_of_spins
    int size_squared_squared = size_squared * size_squared;

    double EE = E * E; double MM = M * M;

    boltzmann(0)  = 	exp(8.0 * beta);	// -8 +8
    boltzmann(4)  = 	exp(4.0 * beta);	// -4 +8
    boltzmann(8)  = 	1;					//  0 +8
    boltzmann(12) = 	exp(-4.0 * beta);	//  4 +8
    boltzmann(16) = 	exp(-8.0 * beta);	//  8 +8

    int cycle;
    for (cycle = 0; cycle <= num_of_mc_cycles; cycle++) {

        metropolis();

        exp_mag_etc_values(0) += E  / size_squared;             // Energy
        exp_mag_etc_values(1) += EE / size_squared_squared;     // Energy squared

        exp_mag_etc_values(2) += M  / size_squared;             // Mag
        exp_mag_etc_values(3) += MM / size_squared_squared;     // Mag squared

        exp_mag_etc_values(4) += fabs(M) / size_squared;        // abs(Mag)

        // This is for choising how many samples you want to be printed in your files
        if ( cycle % 1000 == 0) {
        //if(cycle == 1 || cycle == 100 || cycle == 10000 || cycle == 1000000) {
            cout << cycle << endl;
            if(to_file) {
                get_lattice_to_file(ofile_spin);
                simulate_to_file(cycle + 1, oflie_expect);
                //ofile_except << num_of_mc_cycles << " " << to_file(0) << " " << to_file(1) << " " << to_file(2) << " " << to_file(3) << " " << to_file(4) << endl;
                //ofile_except << num_of_mc_cycles << "," << to_file(0) << "," << to_file(1) << "," << to_file(2) << "," << to_file(3) << "," << to_file(4) << endl; // .csv

            }
        }
    }
}


void IsingModel::reset_expectations() {

    for(int i = 0; i <= 4; i++)
        exp_mag_etc_values(i) = 0;
}

void IsingModel::simulate_to_file(int num_of_mc_cycles, ofstream &ofile_except) {
/*
    //stringstream stream;
    stream << fixed << setprecision(4) << temp;

    //ofile_spin.open("data/" + this->filename + "_" + to_string(stream) + "_spin.txt");
    ofstream ofile_spin("/data/" + this->filename + "_" + to_string(stream) + "_spin.txt");


    //oflie_expect.open("data/" + this->filename + "_" + to_string(stream) + "_expect.txt");
    ofstream oflie_expect("/data/" + this->filename + "_" + to_string(stream) + "_expect.txt");
*/

    // Expectations/Magnetizations printing to file
    vec to_file = exp_mag_etc_values/(num_of_mc_cycles);
    ofile_except << num_of_mc_cycles << " " << to_file(0) << " " << to_file(1) << " " << to_file(2) << " " << to_file(3) << " " << to_file(4) << endl;
    //ofile_except << num_of_mc_cycles << "," << to_file(0) << "," << to_file(1) << "," << to_file(2) << "," << to_file(3) << "," << to_file(4) << endl; // .csv

}

void IsingModel::time_to_file(bool start) {
    /*
    if(start) //double begin = clock();
        high_resolution_clock::time_point begin = high_resolution_clock::now(); // Start clock

    if(!start) {
        //double end  = clock();
        //time = (end - begin) / CLOCKS_PER_SEC;
        high_resolution_clock::time_point end = high_resolution_clock::now(); // Stop clock
        auto time = duration_cast<nanoseconds>(end.count() - begin.count()) / 1e9; // time difference in nano-sec

        cout << time << endlN
        //cout << setw(15) << (end - begin) / CLOCKS_PER_SEC <<endl;
    */

}

