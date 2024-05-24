#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <time.h>

using namespace std;
const int N=100;
const int M=N-2;
const double PI = 3.14159265358979323846264;

int main()
{
    srand(time(NULL));
    for(int i = 0; i<100; i++)
    {
        cout << rand() << " ";
    }
}