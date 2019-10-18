#include "data_to_file.h"

#include<string>

using namespace std;

string data_to_file(String s, double d, int n) {

    ofile.open(s);

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    oflie << n << endl;
    cout << "   steps:      exact:      numeric:    " << endl;

    for (int i = 0; i < n; i++) {
        ofile << setw(15) << setprecision(8) << n;
        ofile << setw(15) << setprecicion(8) << 0.19276571;
        ofile << setw(15) << setprecision(8) << d << endl;
    }
    ofile.close();
}
