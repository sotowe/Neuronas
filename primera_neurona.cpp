//Libraries
#include<iostream>
#include<cstdlib>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

//FUNCTIONS
double Prob(double t);
double dVdt (double t, double V, double tj);
double rungeKutta(double x0, double y0, double h, double tj);


//main function
int main(void)
{
    // declaration of variables
    double tf, t0, V10, V20, t1, t2, Vthres, Vrest, h, t;
    int n,i;
    V10 = -75;
    V20 = -80;
    Vthres = -70;
    Vrest = -80;
    t0 = 0;
    tf = 10;
    h = 0.001;
    n = ((tf-t0)/h);
    t2 = 0;
    t1 = 0;

    ofstream file1;
    file1.open("dosneuronas.txt");

    // We make the iteration, combining both neurons
    for (i=1; i<=n; i++)
    {
        V10 = rungeKutta(t0,V10, h, t2);
        if (V10 >=Vthres)
        {
            V10 = Vrest;
            t1 = t0;
        }

        V20 = rungeKutta(t0, V20, h, t1);
        if (V20 >=Vthres)
        {
            V20 = Vrest;
            t2 = t0;
        }
        t0 = t0 + h;
        file1 << t0 <<"   "<< V10 << "  " << V20 << endl;

    }
file1.close();
return 0;
}


// FUNCTIONS
double Prob(double t)
//function which gives the probability for the channels to be open.
// Input arguments:
//  t: time
// Output arguments:
//  P: probability
{
    // Declaration of variables
    double Pmax, P, ts;

    Pmax = 1;
    ts = 0.5;

    P = 1/ts*exp(1-t/ts);
    return P;
}

double dVdt (double t, double V, double tj)
// Function that calculate the differential equation of the I&F for coupled neurons
{
    double EL, ts, g, dV, tm, Es, Rm, Ie;
    EL = - 60;
    ts = 1;
    g = 0.01;
    tm = 1;
    Es = 0;
    Rm = 1;
    Ie = 1;

    dV = EL - V - g*Prob(t-tj)*(V - Es) - Rm*Ie;

    return dV/tm;
}


double rungeKutta(double x0, double y0, double h, double tj)
// Count number of iterations using step size or
// step height h
{

    double k1, k2, k3, k4, k5;
    double y;

    // Apply Runge Kutta Formulas to find
    // next value of y
    k1 = h*dVdt(x0, y0, tj);
    k2 = h*dVdt(x0 + 0.5*h, y0 + 0.5*k1, tj);
    k3 = h*dVdt(x0 + 0.5*h, y0+ 0.5*k2, tj);
    k4 = h*dVdt(x0 + h, y0 + k3, tj);

    // Update next value of y
    y = y0 + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);;

    return y;
}
