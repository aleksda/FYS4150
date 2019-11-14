#include "isingmodel.hpp"

#include <armadillo>
#include <iomanip>
#include <random>
#include <string>

using namespace std; using namespace arma;

stringstream stream;

/*
ofstream ofile_spin;   ofstream oflie_expect;
ofstream ofile_energy; ofstream ofile_accept;
*/

IsingModel::IsingModel(int size, int data, string filename) {

    this->size     = size;
    this->data     = data;
    this->filename = filename; // + ".csv"

    int num_of_spins = size * size;
    int E, M = 0;

    // vectors & matrices
    spins = Mat<int>(size, size); // Mat<> before spins not working?
    index = Col<int>(size + 2); // FILL WITH ARMA ?

    boltzmann  = vec(17, fill::zeros); // CHANGE to arma vector ?
    exp_value  = vec(5, fill::zeros); // SAME

    // Random number generators
    // From std library
    mt19937 rand_nums(clock()); // CHANGE to high precission clock ?
    uniform_real_distribution<double> uniform_zero_to_one(0.0, 1.0);
    uniform_int_distribution<int> zero_to_L_distribution(0, size-1); // DEL

    //index(0) = size - 1;
    //for(int i = 0; i < size + 2; i++)
        //index(i) = i - 1;

    //index(size + 1) = 0;

    index = Col<int>(size+2);  // DEL
    for (int i = 1; i < size+2; i++) // DEL
    {
      index(i) = i-1;
    }
    index(0) = size-1; // DEL
    index(size+1) = 0; // DEL
}

void IsingModel::initialize_random() {

    //spins.randu()(); // Not working?
    int pm_one = 0;

    for(unsigned int i = 0; i < this->size; i++) {
        for(unsigned int j = 0; j < this->size; j++) {
            if(uniform_zero_to_one(rand_nums) >= 0.5) pm_one = 1;
            else pm_one = -1;
            spins(i, j) = pm_one;
        }
    }
}

void IsingModel::initialize_up() {
    spins.ones();
}

// Not part of project
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

// Metropolis-Hastings Algorithm
double IsingModel::metropolis() {
/*
    int delta_E;
    int rand_x; int rand_y;

    for(int i = 0; i < num_of_spins; i++) {
        rand_x = rand() % (size-1); rand_y = rand() % (size-1);

        delta_E = get_site_energy(rand_y, rand_x);

        if(delta_E < 0 || (double) rand() <= boltzmann(delta_E + 8)) {

            E += delta_E;
            M -= 2 * spins(index(rand_y), index(rand_x));

            spins(index(rand_y), index(rand_x)) *= 1; // Declare it above the loop?
        }

    }*/
    // DEL ALL THIS
    int delta_E;
    int x_i;
    int y_i;
    int neighbour_sum;

    for (int i = 0; i < size*size; i++)
    {
      //choose random spin
      x_i = zero_to_L_distribution(rand_nums)%size+1; //we plus with one since we will use them on "index_vector", which is shifted by +1
      y_i = zero_to_L_distribution(rand_nums)%size+1; //we plus with one since we will use them on "index_vector", which is shifted by +1

      //calculates the sum of the neihbours
      neighbour_sum = spins( index(y_i), index(x_i+1) ) + spins( index(y_i), index(x_i-1) ) + spins( index(y_i+1), index(x_i) ) + spins( index(y_i-1), index(x_i) );

      //uses neighbour_sum to calculate change in energy
      delta_E = 2*neighbour_sum * spins( index(y_i), index(x_i) ); //J=1

      //metropolis check if we want to flip or not
      if (delta_E < 0 || uniform_zero_to_one(rand_nums) <= boltzmann(delta_E + 8))
      {
        E += delta_E; //update energy
        M -= 2*spins( index(y_i), index(x_i));  //update magnetization
        spins( index(y_i) , index(x_i) ) *= -1; //flip spin
      }
    }
  }

//}

void IsingModel::equilibrate(double temp, int num_of_mc_cycles) { // Better name?

    double beta = 1.0 / temp;

    // MAKE THIS IN A FOR LOOP
    boltzmann(0)  = 	exp(8.0 * beta);	// -8 +8
    boltzmann(4)  = 	exp(4.0 * beta);	// -4 +8
    boltzmann(8)  = 	1;					//  0 +8
    boltzmann(12) = 	exp(-4.0 * beta);	//  4 +8
    boltzmann(16) = 	exp(-8.0 * beta);	//  8 +8

    int variable = num_of_mc_cycles / 10; // Maybe finding a more fitting name for this?

    for (int i = 0; i < variable; i++) {
        metropolis();
    }
}

void IsingModel::simulate(double temp, int num_of_mc_cycles, bool to_file) {

    E = find_energy(); M = find_magnetization();

    int denominator = num_of_mc_cycles / this->data;
    double beta = 1.0 / temp;
/*
    if(to_file)
        simulate_to_file(temp);
*/
    //if(to_file) {
        stream << fixed << setprecision(4) << temp;

        //ofile_spin.open("data/" + this->filename + "_" + stream.str() + "_spin.txt");
        ofstream ofile_spin(this->filename + "_" + stream.str() + "_spin.txt");
        if (!ofile_spin.is_open())
           cout << "Could not open spins file" << endl;

        //oflie_expect.open("data/" + this->filename + "_" + stream.str() + "_expect.txt");
        ofstream oflie_expect(this->filename + "_" + stream.str() + "_expect.txt");
        if (!oflie_expect.is_open())
           cout << "Could not open expect file" << endl;

        //ofile_energy.open("data/" + this->filename + "_" + stream.str() + "_energy.txt");
        ofstream ofile_energy(this->filename + "_" + stream.str() + "_energy.txt");
        if (!ofile_energy.is_open())
           cout << "Could not open energy file" << endl;

        get_lattice_to_file(ofile_spin);
    //}

    // Trick for reducing the FLOPs a little during computation
    int size_squared 		 = size * size;	// Alternativly use num_of_spins
    int size_squared_squared = size_squared * size_squared;

    double EE = E * E; double MM = M * M;

    // MAKE THIS IN A FOR LOOP
    boltzmann(0)  = 	exp(8.0 * beta);	// -8 +8
    boltzmann(4)  = 	exp(4.0 * beta);	// -4 +8
    boltzmann(8)  = 	1;					//  0 +8
    boltzmann(12) = 	exp(-4.0 * beta);	//  4 +8
    boltzmann(16) = 	exp(-8.0 * beta);	//  8 +8

    int cycle;
    for (cycle = 0; cycle <= num_of_mc_cycles; cycle++) {

        metropolis();

        exp_value(0) += E  / size_squared;
        exp_value(1) += EE / size_squared_squared;

        exp_value(2) += M  / size_squared;
        exp_value(3) += MM / size_squared_squared;

        exp_value(4) += fabs(M) / size_squared;

        if ( (cycle) % denominator == 0) {
            cout << cycle << endl;
            if(to_file) {
                ofile_energy << E << "\n";
                get_lattice_to_file(ofile_spin);
                //get_exp_to_file(cycle + 1, outfile_expect);
            }
        }
    }
}

void IsingModel::reset_expectations() {
/*
  exp_value(0) = 0; exp_value(1) = 0;
  exp_value(2) = 0; exp_value(3) = 0;
  exp_value(4) = 0;
*/
    for(int i = 0; i <= 4; i++)
        exp_value(i) = 0;
}

void IsingModel::simulate_to_file(double temp) {
/*
    //stringstream stream;
    stream << fixed << setprecision(4) << temp;

    //ofile_spin.open("data/" + this->filename + "_" + stream.str() + "_spin.txt");
    ofstream ofile_spin("/data/" + this->filename + "_" + stream.str() + "_spin.txt");
    if (!ofile_spin.is_open())
       cout << "Could not open spins file" << endl;

    //oflie_expect.open("data/" + this->filename + "_" + stream.str() + "_expect.txt");
    ofstream oflie_expect("/data/" + this->filename + "_" + stream.str() + "_expect.txt");
    if (!oflie_expect.is_open())
       cout << "Could not open expect file" << endl;

    //ofile_energy.open("data/" + this->filename + "_" + stream.str() + "_energy.txt");
    ofstream ofile_energy("/data/" + this->filename + "_" + stream.str() + "_energy.txt");
    if (!ofile_energy.is_open())
       cout << "Could not open energy file" << endl;
*/

/*
    ofile_accept.open("data/" + this->filename + "_" + stream.str() + "_accepts.txt");
    if (!ofile_accept.is_open())
       cout << "Could not open accept file" << endl;
*/
}
