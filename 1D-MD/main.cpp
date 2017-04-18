#include "main.h"


void calculate_force_component(double& Fx, double& mass, double& omega, double& X)
{
    Fx=-mass*pow(omega,2)*X;
}

void calculate_energy(double& T, double& V, double& mass, double& omega, double& X, double& Vx)
{
    T=0.5*mass*pow(Vx,2);
    V=0.5*mass*pow(omega,2)*pow(X,2);
}

void calculate_position(double& X, double& Xo,double& Fx, double& mass, double& Vx, double& omega, double& time, double& dt)
{
   /*
   if (time == 0)
   {
       Vx=10.0;
       Xo=X-Vx*dt;
   }
   else
   {
      // Do Nothing
   }
   */
   calculate_force_component(Fx, mass, omega, X);
   X=2*X-Xo+Fx/mass*pow(dt,2);
   Vx=(X-Xo)/(2.0*dt);
}

void velocity_verlet(double& X, double& Vx, double& Fx, double& mass, double& omega,  double& dt)
{
    double xnew=0.0, vnew=0.0, fxnew;
    xnew=X+dt*Vx+pow(dt,2)/(2.0*mass)*Fx;
    calculate_force_component(fxnew, mass, omega, xnew);
    vnew=Vx+dt/(2.0*mass)*(Fx+fxnew);
    X=xnew;
    Vx=vnew;
    Fx=fxnew;
}

void md_step(double& T, double& V, double& X, double& Fx, double& mass, double& Vx, double& dt, double& omega)
{
    velocity_verlet(X, Vx, Fx, mass, omega, dt);
    calculate_energy(T,V,mass,omega,X,Vx);
}

int main(int argc, char *argv[])
{
    double mass=1.0, omega=.2;
    double X=1.0, Vx=0.0, Fx=0.0;
    double time = 0.0, dt = 1.0;
    double T = 0.0, V = 0.0, E=0.0;
    int nsteps=10000;

    dt=stod(argv[1]);
    ofstream logfile;
    logfile.open("log.harm");
    calculate_force_component(Fx, mass, omega, X);
    calculate_energy(T, V, mass, omega, X, Vx);
    E=T+V;
    //logfile << "TIME X VX FX T V E" << endl;
    logfile << time << " " << X << " " << Vx << " " << Fx << " " << T << " " << V << " " << E <<endl;
    for (int t=1; t<nsteps; t++)
    {
        time = t*dt;
        md_step(T, V, X, Fx, mass, Vx, dt, omega);
        E=T+V;
        logfile << time << " " << X << " " << Vx << " " << Fx << " " << T << " " << V << " " << E <<endl;
    }
    

    return 0;
}
