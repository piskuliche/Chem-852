#include "main.h"

//This is the main program file for the third computational exercise  (problem 2) of Chem-852.
//All Harmonic Oscillator parameters are taken from McQuarrie Quantum Mechanics and NIST
//All Morse Potential parameters are taken from
//  D. Konowalow, J. Herschfelder, "Morse Potential Parameters for O-O, N-N, and N-O Interactions", Phys. Fluids, 4, (5), 637-642, (1960)


//This calculation assumes a room temperature of 298 K.


const double K=1142;                  // N m^-1   
const double Req=0.0;             // Angstroms
const double Beta=2.78;               // Angstroms^-1
const double A = sqrt(48);            // kcal/(mol * Angstrom^2)
const double B = 1;                   // kcal/(mol * Angstrom^4)
const double De=5.211;                // Ev
const double JtoEv=6.242*pow(10,18);  // Ev
const double AngtoM=1*pow(10,-10);    // m
const double BoltzBeta=38.9413;       // Ev^-1
const double global_minimum=-5;    // Minimum for histogram
const double bin_separation=0.01;      // How to count histogram bins
const double GAUSSHEIGHT=-12.0;       //Height of Umbrella potential  
const double GAUSSPOS=0.0;            //Position of Umbrella Potential
const double GAUSSWIDTH=0.3;          //Width of Umbrella Potential
int shift_bin = abs(global_minimum/bin_separation);
const double kcaltoev = .043363;
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
double umbrella_pot(double& R)
{
    double RSEP=calc_sep(R);
    double umbV = kcaltoev*GAUSSHEIGHT*exp(-pow((RSEP-GAUSSPOS),2)*GAUSSWIDTH);
    return umbV;
}

double calc_energy(int& en_type, double& R)
{
    double current_energy = 0;
    double r_sep=calc_sep(R);
    if (en_type == 0) //Harmonic Oscillator Parameters
    {
        current_energy = 0.5 * K * pow(r_sep,2) * JtoEv * pow(AngtoM,2);
    }
    else if (en_type == 1)
    {
        current_energy = -0.5 * A * pow(r_sep,2) + 0.25 * B * pow(r_sep,4) + pow(A,2)/(4*B);
        current_energy = current_energy*kcaltoev;
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

void histogram_data(bool& umbsampleflag, int print, double& value, vector<int>& histogram_storage, int& MC_STEPS)
{
    if (print == 0)
    {
        //Calculates a bin number
        int bin_no = (value - global_minimum)/(bin_separation);
        //Adds a count to the bin number mentioned.
        histogram_storage[bin_no] += 1;
    }
    else
    {
        ofstream hist;
        hist.open("metropolis.hist");
        int num_bins = static_cast<int>(2*shift_bin);
        for (int i=0; i<num_bins; i++)
        {
            if (umbsampleflag==0)
            {
                hist << global_minimum+i*bin_separation << " " << histogram_storage[i] << endl;
            }
            else
            {
                double potpoint=global_minimum+i*bin_separation;
                int count =0;
                for (int m=0; m<histogram_storage.size(); m++)
                        {
                            count += histogram_storage[m];
                        }

                double weightbin=histogram_storage[i]/(-umbrella_pot(potpoint)*count);
                hist << global_minimum+i*bin_separation << " " << histogram_storage[i] << " " << weightbin << endl;
            }
        }
        hist.close();
    }
}


void metropolis_mc_step(bool umbsampleflag, vector<int>& histogram_storage, double& r_old, double& r_new, double& dr, double& U_Old, double& U_New, int& en_type, int& acc, int& att, double& mc_r, double& mc_r2, double& mc_energy, int& MC_STEPS)
{
    //Call Shift algorithm to cause a displacement in r
    shift(r_old, r_new, dr);
    //Calculate the New Energy
    if (umbsampleflag == 0)
    {
        U_New = calc_energy(en_type,r_new);
    }
    else
    {
        U_New = calc_energy(en_type, r_new) + umbrella_pot(r_new);
    }
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
    
    //Addition for the histogram subroutine
    double tmpsep=r_old-Req;
    int tmpprnt = 0;
    histogram_data(umbsampleflag, tmpprnt, tmpsep, histogram_storage, MC_STEPS);
}


int main()
{
    srand(time(NULL));
    //Parameter Variables
    double dr = 0;
    int bins = static_cast<int>(2*shift_bin);
    vector<int> histogram_storage(bins,0);
    int en_type = 0;
    bool umbsampleflag = 0; 
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
    cout << "Please enter 0 for harmonic, 1 for a quartic potential, or any other key for the morse potential" << endl;
    cin >> en_type;
    cout << "Please enter a number of steps to run" << endl;
    cin >> MC_STEPS;
    cout << "Please enter a non-zero dumpfrequency" << endl;
    cin >> dumpfreq; 
    cout << "Would you like to use the umbrella sampling technique? (0 for no, 1 for yes)" << endl;
    cin >> umbsampleflag;
    cout << "Running" << endl;
    if (umbsampleflag == 1)
    {
        cout << "Beginning Umbrella Sampling" << endl;
    }
    
    ofstream potential;
    potential.open("harm.log");
    double dist=-0.5;
    double px =0;
    while (dist < 0.5)
    {
        px = sqrt(2*3.14159)/sqrt(K*BoltzBeta*AngtoM*JtoEv)*exp(-BoltzBeta*calc_energy(en_type,dist));
        potential << dist << " " << px << endl;
        dist += .001;
    }



    ofstream log;
    log.open("mc_1d.log");
    ofstream output;
    output.open("mc_1d.out");
    log << "Run_Parameters:" << endl;
    log << "dr = " << dr << endl;
    log << "mc_steps = " << MC_STEPS << endl;
    output << "STEP ACCRATIO <dR> <dR2> <ENERGY>" << endl;
    //Initialize r
    /*
    if (random_generator() < 0.5)
    {
        r_old=Req + random_generator()*Req;
    }
    else
    {
        r_old=Req - random_generator()*Req;
    }
    */
    //Init in minimum
    r_old=2.6;
    if (umbsampleflag == 0)
    {
        U_Old = calc_energy(en_type,r_old);
    }
    else
    {
        U_Old = calc_energy(en_type, r_old) + umbrella_pot(r_old);
    }
    for (int i = 0; i < MC_STEPS; i++)
    {
        metropolis_mc_step(umbsampleflag, histogram_storage, r_old, r_new, dr, U_Old, U_New, en_type, acc, att, mc_r, mc_r2, mc_energy, MC_STEPS);
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
    int print=1;
    double tmpval=0.0;
    histogram_data(umbsampleflag, print, tmpval, histogram_storage, MC_STEPS);
    return 0;
}
