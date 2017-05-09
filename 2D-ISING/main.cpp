#include "main.h"
void build_lattice(vector< vector<int> >& GRID, int& rows, int& cols)
//Initializes the Grid Based on rows, cols.
{
    GRID.resize(rows+2);
    for (int i =0; i<=rows+1; i++)
    {
        GRID[i].resize(cols+2,0);
    }
}

void output_grid(vector< vector<int> >& GRID, int& rows, int& cols)
//General Subroutine for Displaying the Grid
{
    cout << "--------------------------" << endl;
    for (int i = 0; i<=rows+1; i++)
        {
            if (i == 1)
            {
                cout << "-  | ";
                for (int j = 1; j <= cols; j++)
                {
                    cout << "- ";
                }
                cout << " | -";
                cout << endl;
            }
            for (int j = 0; j<=cols+1; j++)
            {
                if (j == 1)
                {
                    cout << " | ";
                }
                cout << GRID[i][j] << " ";
                if (j == cols)
                {
                    cout << " | ";
                }
            }
            cout << endl;
            if (i==rows)
            {
                cout << "-  | ";
                for (int j = 1; j <= cols; j++)
                {
                    cout << "- ";
                }
                cout << " | -";
                cout << endl;
            }
        }
    cout << "--------------------------" << endl;
}

double random_generator()
//Calculates a random number on the range[0:1]
{
    double random_number;
    random_number = rand() / (double) RAND_MAX;
    return random_number;
}

void insertions(vector< vector<int> >& GRID, int& rows, int& cols)
//Subroutine which inserts on each lattice point based on the probability
{
    for (int i = 1; i<=rows; i++)
    {
        for (int j = 1; j <= cols; j++)
        {
            if (random_generator() < 0.5)
            {
                GRID[i][j]=3;
            }
            else
            {
                GRID[i][j]=1;
            }
        }
    }
    for (int j = 1; j <= cols; j++)
    {
        GRID[0][j]=GRID[rows][j];
        GRID[rows+1][j]=GRID[1][j];
    }
    for (int i = 1; i <= rows; i++)
    {
        GRID[i][0]=GRID[i][cols];
        GRID[i][cols+1]=GRID[i][1];
    }
}

double energy_calc(vector< vector<int> >& lattice, int& rows, int& cols, int& i, int&j)
{
    double ENERGY=0.0;
    ENERGY=lattice[i-1][j]+lattice[i+1][j]+lattice[i][j-1]+lattice[i][j+1]-2*4;
    return ENERGY;
}
double calculate_en_change(vector< vector<int> >& lattice, int& rows, int& cols, int& i, int& j, double& J)
{
    double DELE=0.0;
    if ( lattice[i][j] == 1)
    {
        DELE=-2*J*energy_calc(lattice,rows,cols,i,j);
    }
    else
    {
        DELE=2*J*energy_calc(lattice,rows,cols,i,j);
    }
    return DELE;
}

int select_random_index(int& dim)
{
    double r1=0.0;
    r1=random_generator();
    r1=r1*dim;
    int choose_out=(int) ceil(r1);
    return choose_out;
}

void metropolis_step(vector< vector<int> >& lattice, int& rows, int& cols, double& EOLD, int& acc, int& att, double& J, double& beta)
{
    int cycstep=rows*cols;
    for (int nstep=0; nstep < cycstep; nstep++)
    {
        int tmp_i=select_random_index(rows);
        int tmp_j=select_random_index(cols);
        double DELE=calculate_en_change(lattice, rows, cols, tmp_i, tmp_j, J);
        double tmpE=EOLD+DELE;
        if (DELE < 0)
        {
            //Accept
            EOLD=tmpE;
            if (lattice[tmp_i][tmp_j] == 1)
            {
                lattice[tmp_i][tmp_j]+=2;
            }
            else
            {
                lattice[tmp_i][tmp_j]-=2;
            }
            acc+=1;
            att+=1;
        }
        //fix this
        else if (exp(-DELE/beta) > random_generator())
        {
            //Secondary Accept
            EOLD=tmpE;
            if (lattice[tmp_i][tmp_j] == 1)
            {
                lattice[tmp_i][tmp_j]+=2;
            }
            else
            {
                lattice[tmp_i][tmp_j]-=2;
            }
            acc+=1;
            att+=1;
        }
        else
        {
            // Reject
            att+=1;
        }
        
       if (tmp_i == 1)
       {
           lattice[rows+1][tmp_j]=lattice[tmp_i][tmp_j];
       }
       else if (tmp_i == rows)
       {
           lattice[0][tmp_j]=lattice[tmp_i][tmp_j];
       }
       if (tmp_j == 1)
       {
           lattice[tmp_i][cols+1]=lattice[tmp_i][tmp_j];
       }
       else if (tmp_j == cols)
       {
           lattice[tmp_i][0]=lattice[tmp_i][tmp_j];
       }
    }
    //END LOOP
    //CALL FULL E CALC
}

double calculate_av_spin(vector< vector<int> >& lattice, int& rows, int& cols)
{
    int N=rows*cols;
    double SUM=0.0;
    for (int i=1; i <= rows; i++)
    {
        for (int j = 1; j<=cols; j++)
        {
            SUM+=(lattice[i][j]-2.0);
        }
    }
    double av=0.0;
    av=SUM/(double) N;
    return av;
}

int calc_dist(int& rows, int& cols, int&i, int& j, int& i_in, int& j_in)
{
    double del_i=i_in-i, del_j=j_in-j;
    double Lr=rows/2.0, Lc=cols/2.0;
    if (del_i > Lr)
    {
        del_i-=Lr;
    }
    else if (del_i < -Lr)
    {
        del_i+=Lr;
    }
    if (del_j > Lc)
    {
        del_j-=Lc;
    }
    else if (del_j < -Lc)
    {
        del_j+=Lc;
    }
    double dist=sqrt(pow(del_i,2)+pow(del_j,2));
    int dist_bin=floor(dist/(Lr/10.0));
    return dist_bin;
}


void calculate_av_corr(vector< vector<int> >& lattice, vector<int>& n_r, vector<double>& c_r, vector<double>& a_r,vector<double>& ss_r, int& rows, int& cols, double& avspin)
{
    double SUM=0.0;
    int bin=0;
    for (int i=1; i<=rows-1; i++)
    {
        for (int j=1; j<=cols-1; j++)
        {
            for (int i_in=i+1; i_in<=rows; i_in++)
            {
                for (int j_in=j+1; j_in<=cols; j_in++)
                {
                    bin=calc_dist(rows, cols, i, j, i_in, j_in);
                    n_r[bin]+=1;
                    c_r[bin]+=(lattice[i_in][j_in]-2)*(lattice[i][j]-2);
                    a_r[bin]=c_r[bin]/(double) n_r[bin];
                }
            }
        }
    }
    for (int r=0; r<10; r++)
    {
        ss_r[r]=a_r[r]-pow(avspin,2);
    }
}

int main(int argc,char *argv[])
{
    int rows = 20, cols = 20, ncycl=1000;
    double T=stod(argv[1]);
    int acc = 0, att = 0;
    double kb=1.0;//0.0019872041;
    double beta=kb*T;
    double EOLD=0.0;
    double Tc=500.0, J=0.0;
    double avspin=0.0;
    //J=Tc*kb/4.0;
    J=1.0;
    srand(time(NULL));
    vector< vector<int> > lattice;
    vector<double> c_r(100,0), a_r(100,0), ss_r(100,0);
    vector<int> n_r(100,0);
    build_lattice(lattice,rows,cols);
    output_grid(lattice,rows,cols);
    insertions(lattice,rows,cols);
    output_grid(lattice,rows,cols);
    string logfile="log.ising."+to_string(T)+".out";
    string corfile="ss_r."+to_string(T)+".out";
    ofstream output;
    output.open(logfile);
    for (int cyc=0; cyc < ncycl; cyc++)
    {
        metropolis_step(lattice, rows, cols, EOLD, acc, att, J, beta);
        avspin=calculate_av_spin(lattice, rows, cols);
        calculate_av_corr(lattice, n_r, c_r, a_r, ss_r, rows, cols, avspin);
        output << cyc << " "  << avspin << endl;
        cout << avspin << " " <<att << " " << acc << " " << (float) (acc/(float)att) << endl;
        output_grid(lattice, rows, cols);
    }
    output.close();
    ofstream corrout;
    corrout.open(corfile);
    for (int r=0; r<10; r++)
    {
        corrout << r << " " << ss_r[r]<< " "  << a_r[r] << " " << n_r[r] << " " << c_r[r] << endl;
    }
    //OUTPUT STUFF
    return 0;
}
