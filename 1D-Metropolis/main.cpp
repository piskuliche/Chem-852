#include "main.h"

//This is the main program file for the second computational exercise of Chem-852.
//All Harmonic Oscillator parameters are taken from McQuarrie Quantum Mechanics and NIST
//All Morse Potential parameters are taken from
//  D. Konowalow, J. Herschfelder, "Morse Potential Parameters for O-O, N-N, and N-O Interactions", Phys. Fluids, 4, (5), 637-642, (1960)


//This calculation assumes a room temperature of 298 K.


const double K=1142;                  // N m^-1   
const double Req=1.20752;             // Angstroms
const double Beta=2.78;               // Angstroms^-1
const double De=5.211;                // Ev
const double JtoEv=6.242*pow(10,18);  // Ev
const double AngtoM=1*pow(10,-10);    // m
const double BoltzBeta=38.9413;       // Ev^-1

double random_generator()
{
    double random_number;
    random_number = rand() / (double) RAND_MAX;
    return random_number;
}
double calc_sep(double& R)
{
    double separation=0;
    separation=R-Req;
    return separation;
}
double calc_energy(bool& en_type, double& R)
{
    double current_energy = 0;
    double r_sep=calc_sep(R);
    if (en_type == 0) //Harmonic Oscillator Parameters
    {
        current_energy = 0.5 * K * pow(r_sep,2) * JtoEv * pow(AngtoM,2);
    }
    else
    {
        current_energy = De*pow((1-exp(-Beta*(r_sep))),2);
    }
    return current_energy;
}
void shift(double& r_old, double& r_new,double& dr)
{
    if (random_generator()<0.5)
    {
        r_new = r_old + dr*random_generator();
    }
    else
    {
        r_new = r_old - dr*random_generator();
    }
}

void metropolis_mc_step(double& r_old, double& r_new, double& dr, double& U_Old, double& U_New, bool& en_type, int& acc, int& att, double& mc_r, double& mc_r2, double& mc_energy)
{
    //Call Shift algorithm to cause a displacement in r
    shift(r_old, r_new, dr);
    //Calculate the New Energy
    U_New = calc_energy(en_type,r_new);
    //Test for Acceptance and Rejectance Criteria
    if (U_New < U_Old) //Accept
    {
        r_old=r_new;
        U_Old=U_New;
        acc += 1;
        att += 1;
    }
    else
    {
        double Boltzmann_Factor=exp(-BoltzBeta*(U_New-U_Old));
        double rand_val=random_generator();
        if (Boltzmann_Factor > rand_val) //Accept
        {
            r_old=r_new;
            U_Old=U_New;
            acc += 1;
            att += 1;
        }
        else //Reject
        {
            att += 1;
        }
    }
    //This section ensures that counting is taken care of at each step.
    mc_r        += (r_old-Req);
    mc_r2       += pow((r_old-Req),2);
    mc_energy   += U_Old;
}

int main()
{
    srand(time(NULL));
    //Parameter Variables
    double dr = 0;
    bool en_type = 0;
    int dumpfreq = 0;
    int MC_STEPS = 0;
    //Progress Variables
    double U_Old=0;
    double U_New=0;
    double r_old=0;
    double r_new=0;
    //Output variables
    int acc=0;
    int att=0;
    double mc_r=0;
    double mc_r2=0;
    double mc_energy=0;

    cout << "Welcome to Metropolis Monte Carlo in 1D" << endl;
    cout << "Please enter a value of dr" << endl;
    cin >> dr;
    cout << "Please enter 0 for harmonic or 1 for morse potential" << endl;
    cin >> en_type;
    cout << "Please enter a number of steps to run" << endl;
    cin >> MC_STEPS;
    cout << "Please enter a non-zero dumpfrequency" << endl;
    cin >> dumpfreq; 
    cout << "Running" << endl;

    ofstream log;
    log.open("mc_1d.log");
    ofstream output;
    output.open("mc_1d.out");
    log << "Run_Parameters:" << endl;
    log << "dr = " << dr << endl;
    log << "mc_steps = " << MC_STEPS << endl;
    output << "STEP ACCRATIO <dR> <dR2> <ENERGY>" << endl;
    //Initialize r
    if (random_generator() < 0.5)
    {
        r_old=Req + random_generator()*Req;
    }
    else
    {
        r_old=Req - random_generator()*Req;
    }
    U_Old = calc_energy(en_type, r_old);
    for (int i = 0; i < MC_STEPS; i++)
    {
        metropolis_mc_step(r_old, r_new, dr, U_Old, U_New, en_type, acc, att, mc_r, mc_r2, mc_energy);
        if (i%dumpfreq == 0)
        {
            cout << "Beginning step " << i << endl;
            output << i << " " << acc / ((float) att) << " " << mc_r / ((float) att) << " " << mc_r2 / ((float) att) << " " << mc_energy / ((float) att) << endl;
        }
    }
    log << "Final dr = " << dr << endl;
    log << "acceptance ratio = " << acc/((float) att) << endl;
    log << "<dr> = " << mc_r / ((float) att) << endl;
    log << "<dr^2> = " << mc_r2 / ((float) att) << endl;
    log << "<E> = " << mc_energy / ((float) att) << endl;
    log.close();
    cout << "Run Complete" << endl;
    return 0;
}
