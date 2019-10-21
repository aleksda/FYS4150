#include "g_legendre.h"
#include "g_laguerre.h"

#include <iostream>


using namespace std;

int main()
{
    double pi = 3.1415;
    gauss_laguerre(0,2,0,pi,0,2*pi,20);

    return 0;
}
