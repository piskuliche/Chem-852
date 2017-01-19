#include "main.h"

double random_generator()
{
    double random_number;
    random_number = rand() / (double) RAND_MAX;
    return random_number;
}




int main()
{
    for (int i=0; i<10; i++)
    {
        double r = random_generator();
        cout << r << endl;
    }
    return 0;
}
