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

int main(void)
{
    int N; //Number of neurons
    int i, j, k, counterr, counterV; // counters
    int sum; // Number of neurons for the average of the potential.
    int connectivity [10000]; //Number of conexions of every Neuron;
    int neuron_connected; //number of the connected neuron
    double aux, h, /*time increment*/ t, /*time*/tmax, tmin, s;//mean synaptic activation
    double Vp, /*peak*/ tau, average, averageV, averager, /*average potential*/ rate; /*firing rate*/
    double Vthres, /*Vthreshold*/ Vr; /*Vrest*/
    double etamedia, dt, /*time we use for rate*/ Ndt; //N*dt
    double J; //Synaptic weight
    double Js; //J*s
    double taverage;
    double prob_connection; //probability for two neurons to be connected
    double tpot[10000]; // time when the potential was sent
    double eta[10000];
    double refrac_period[10000]; //Refractory period for every neuron.
    double I[10000]; //Input current
    double V[10000]; //Neuron potential
    double rast[300]; //Raster vector of 300 neurons;
    const double pi = atan(1)*4; //PI
    bool spike[10000]; //spiking neurons
    vector< vector<int> > connection; //Vector of connections

    //Generador de aleatorios
    random_device rd;
    mt19937 gen(45454536345331);
    uniform_real_distribution<> con(0,100);

    cout<<"eta"<<"  "<<"   V  "<<"  "<<"   r"<<endl;
    N = 10000;
    J = 20.0;
    tau = 1.0e-3;
    h = 1.0e-4;
    tmax = 15;
    tmin = -15;
    taverage = 10;
    Vp = 100.0;
    Vthres = Vp;
    Vr = -Vp;
    rate = 0;
    dt = 1e-2;
    Ndt = N*dt;
    prob_connection = 0.1;

    ofstream file1;
    ofstream file2;
    ofstream file3;

    file1.open("archivostxt/meanaverageV_J20_ro0_Vo-3_co10-2.txt");
    file2.open("archivostxt/meanrate_J20_ro0_Vo-3_co10-2.txt");
    file3.open("archivostxt/ratecheck.txt");

    //Creamos el número de filas del vector
    for (i = 0; i<N; i++)
    {
        connection.push_back(vector<int>());
    }


    //Creamos las distintas conectividades
    for (i=0;i<N;i++)
        {
            for(j=i+1;j<N;j++)
            {
               if (con(gen)<prob_connection)
               {
                    connection[i].push_back(j);
                    connection[j].push_back(i);
               }

            }
            connectivity[i]= connection[i].size();
        }


    // We generate a provisional refractory period for every Neuron
    for (j=0;j<N;j++)
    {
        refrac_period[j]= 2.0/Vp;
    }

    for (etamedia = -7; etamedia <=-5.8; etamedia+=0.2)
    {
        counterr = 0;
        counterV = 0;
        averageV = 0;
        averager = 0;

        // We calculate eta for each Neuron. Only needed once.
        for (j =0; j < N; j++)
        {
            aux = (pi/2.0)*((2.0*(j+1)-N-1)/(N+1.0));
            eta[j] = etamedia + tan(aux);
        }

        // V are all 0
        for (j=0;j<N; j++)
        {
            V[j]= -3.0;
        }

        // No Neuron has been fired yet
        for (j=0;j<N; j++)
        {
            tpot[j]=-10000;
        }

        rate = 0;
        for (t=tmin; t <= tmax; t += h)
        {

            // Next we calculate the potential of every Neuron using Gauss method for differential equations. Also
            // we calculate the average potential.


            // We calculate the rate each 10^-2 s
            int ftime = t*10000;
            average = 0;
            sum = 0;
            for (j=0;j <N; j++)
            {
            spike[j] = false;

                s = 0;
                //Only if the neuron is not in refractory period.
                if (t - tpot[j] > refrac_period[j])
                {
                    if (t<0)
                    {
                        s = 0;
                    }
                    else
                    {
                        // We calculate s for Neuron j
                        for (i = 0; i < connectivity[j]; i++)
                        {
                            neuron_connected = connection[j][i];
                            if((t - tpot[neuron_connected] <= tau)&&(t - tpot[neuron_connected] >= 0))
                            {
                                s += 1;
                            }
                        }
                        s = 1.0*s/(connectivity[j]*1.0*tau);
                    }
                    Js = J*s;
                    I[j] = eta[j] + Js;
                    Euler(V[j],I[j],h);

                     // If the potential of the neuron reach Vthres, the neuron sends an impulse and goes to a refractory state.
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
                counterV += 1;
            }
            if (ftime%100 == 0)
            {
                rate = 1.0*rate/(Ndt);
                if (t > taverage)
                {
                    averager +=rate;
                    counterr += 1;
                }
                file3 << t <<"  "<< rate << endl;
                rate = 0;
            }

        }

        file1<<etamedia<<"  "<<averageV/counterV<<endl;
        file2<<etamedia<<"  "<<averager/counterr<<endl;
        cout<<etamedia<<"  "<<averageV/counterV<<"  "<<averager/counterr<<endl;
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
