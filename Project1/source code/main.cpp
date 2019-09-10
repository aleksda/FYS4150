//# include "project1b.h"
//# include "project1c.h"
# include "project1d.h"
//# include "project1e.h"
# include <armadillo>
# include <iostream>
# include <string>

using namespace std;
//using namespace arma;

int main(int argc, char* argv[]){
    //int exponentb;
    int exponent;
    string filename;


    // We read also the basic name for the output file and the highest power of 10^n we want
    if( argc <= 1 ){
        cout << "Test " << "\n";
        cout << "Bad Usage: " << argv[0] << " read also file name on same line and max power 10^n for b first and c,d,e second" << endl;
        exit(1);
    }
    else{

        filename = argv[1]; // first command line argument after name of program
        //exponentb = atoi(argv[2]);
        exponent = atoi(argv[2]);

    }
    //project1b(fn, exponent);
    //project1c(filename, exponent);
    project1d(filename, exponent);
    //project1e(filename, exponent);

    return 0;
}
