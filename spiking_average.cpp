//Libraries
#include<iostream>
#include<cstdlib>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;

// FUNCTIONS
double dVdt(double V, double I);
void Euler(double &V, double I, double h);
double errorEstandar(double media2, double media, double medidas);

int main(void)
{
    int N; //Number of neurons
    int i, j, k, counterr, counterV; // counters
    int sum; // Number of neurons for the average of the potential.
    double aux, h, /*time increment*/ t, /*time*/tmax, tmin, s;//mean synaptic activation
    double Vp, /*peak*/ tau, average, averageV, averager, /*average potential*/ rate; /*firing rate*/
    double averageV2, averager2; //sum(V_i^2) and sum(r_i^2)
    double Vthres, /*Vthreshold*/ Vr; /*Vrest*/
    double errorr, errorV;
    double taverage;
    double J; //Synaptic weight
    double Js; //J*s
    double etamedia, dt, /*time we use for rate*/ Ndt; //N*dt
    double tpot[10001]; // time when the potential was sent
    double eta[10001];
    double rast[300]; //Raster vector of 300 neurons;
    double I[10001]; //Input current
    double V[10001]; //Neuron potential
    double refrac_period[10001]; //Refractory period for every neuron.
    const double pi = atan(1)*4; //PI
    bool spike[10001]; //spiking neurons


    cout<<"eta"<<"  "<<"   V  "<<"  "<<"   r"<<"    "<<"errorV"<<"   "<<"errorr"<<endl;
    N = 10000;
    J = 20.0;
    tau = 1.0e-3;
    h = 1.0e-4;
    tmax = 45;
    tmin = -15;
    taverage = 15;
    Vp = 100.0;
    Vthres = Vp;
    Vr = -Vp;
    rate = 0;
    dt = 1e-2;
    Ndt = N*dt;
    ofstream file1;
    ofstream file2;
    ofstream file3;

    file1.open("archivostxt/meanaverageV_J20_ro0_Vo-3_long_2.txt");
    file2.open("archivostxt/meanrate_J20_ro0_Vo-3_long_2.txt");
    file3.open("archivostxt/ratecheck.txt");

    for (j=1;j<=N;j++)
    {
        refrac_period[j]= 2.0/Vp;
    }
    for (etamedia = -4; etamedia <-2.8; etamedia+=0.2)
    {
        counterr = 0;
        counterV = 0;
        averageV = 0;
        averager = 0;
        averageV2 = 0;
        averager2 = 0;

        // We calculate eta for each Neuron. Only needed once.
        for (j =1; j <= N; j++)
        {
            aux = (pi/2.0)*((2.0*j-N-1)/(N+1.0));
            eta[j] = etamedia + tan(aux);
        }

        // V are all 0
        for (j=1;j<=N; j++)
        {
            V[j]= -3;
        }

        // No Neuron has been fired yet
        for (j=1;j<=N; j++)
        {
            tpot[j]=-10000;
        }

        rate = 0;
        for (t=tmin; t <= tmax; t += h)
        {
         // We calculate  the mean synaptic activation for every Neuron. The first 5 seconds we have r = 1;
                if (t<0)
                {
                    s = 0;
                }
                else
                {
                    s = 0;
                    for (j=1; j<=N; j++)
                    {
                        if((t - tpot[j] <= tau)&&(t - tpot[j] >= 0))
                        {
                            s += 1;
                        }

                    }
                    s = 1.0*s/(N*1.0*tau);
                }
                Js = J*s;


            // Next we calculate the potential of every Neuron using Gauss method for differential equations. Also
            // we calculate the average potential.
            average = 0;
            sum = 0;

            // We calculate the rate each 10^-2 s
            int ftime = t*10000;

            for (j=1;j <=N; j++)
            {
            spike[j] = false;
                //Only if the neuron is not in refractory period.
                if (t - tpot[j] > refrac_period[j])
                {

                    I[j] = eta[j] + Js;
                    Euler(V[j],I[j],h);

                     // If the potential of the neuron reach Vthres, the neuron sends an impulse and
                     //goes to a refractory state.
                    if (V[j]>= Vthres)
                    {
                        refrac_period[j]=2.0/V[j];
                        tpot[j] = t + 1.0/V[j];
                        V[j] = -V[j];
                        rate = rate += 1;
                        spike[j] = true;
                    }

                }
                // We calculate the average potential in the network. Only no refractory.
                if (t - tpot[j] > refrac_period[j])
                {
                    average += V[j];
                    sum += 1;

                }
            }
            average = 1.0*average/sum;

            if (t > taverage)
            {
                averageV += average;
                averageV2 += average*average;
                counterV += 1;
            }
            if (ftime%100 == 0)
            {
                rate = 1.0*rate/(Ndt);
                if (t > taverage)
                {
                    averager +=rate;
                    averager2 += rate*rate;
                    counterr += 1;
                }
                file3 << t <<"  "<< rate << endl;
                rate = 0;

            }

        }
        averageV = averageV/counterV;
        averager = averager/counterr;
        errorV = errorEstandar(averageV2, averageV, counterV);
        errorr = errorEstandar(averager2, averager, counterr);
        file1<<etamedia<<"  "<<averageV<<"  "<<errorV<<endl;
        file2<<etamedia<<"  "<<averager<<"  "<<errorr<<endl;
        cout<<etamedia<<"  "<<averageV<<"  "<<averager<<"  "<<errorV<<"  "<<errorr<<endl;
    }
    file1.close();
    file2.close();
    file3.close();
    return 0;

}



// FUNCTIONS

double dVdt(double V, double I)
// Differential equation for QIF neurons
// Input:
//  V: potential
//  Ij: input current. Ij = nj + Js(t) + I'(t)
// Output:
//  dV
{
    double dV;
    dV = V*V + I;
    return dV;
}


void Euler(double &V, double I, double h)
    {
        V = V + h*dVdt(V,I);
        return;
    }

double errorEstandar (double media2, double media, double medidas)
{
    double error;
    error = 1.96*sqrt(1.0*(media2 - medidas*media*media))/medidas;
    return error;
}
