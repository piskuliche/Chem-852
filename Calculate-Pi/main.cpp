#include "main.h"

double random_generator()
{
    double random_number;
    random_number = rand() / (double) RAND_MAX;
    return random_number;
}




int main()
{
    //Initialize random number generator
    srand(time(NULL));
    //Initialize Variables
    double x=0.0, y=0.0, rsq=0.0, pi_calc=0.0, err =0.0;
    double mc_steps = 1*pow(10,12);
    double N_tot=0, N_in=0;
    const double PI=M_PI;
    ofstream output;
    output.open("calcpi.out");
    output.precision(15);
    for (double step=0; step <mc_steps; step+=1.0)
    {
        //Choose an X and Y coordinate.
        x=random_generator(), y=random_generator();
        rsq=pow(x,2)+pow(y,2);
        if (rsq < 1.0)
        {
            N_in+=1;
        }
        N_tot+=1;
        pi_calc = 4 * N_in / N_tot;
        err = (PI-pi_calc)/ PI; 
        if (fmod(step,pow(10,6)) == 0)
        {
            output << "Step " << step+1 << " " << pi_calc << " " << err << endl; 
        }
    }
    
    
    
    return 0;
}
