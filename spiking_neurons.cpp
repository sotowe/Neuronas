//Libraries
#include<iostream>
#include<cstdlib>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include<vector>

using namespace std;


// FUNCTIONS
double dVdt(double V, double I);
void Euler(double &V, double I, double h);

int main(void)
{
    int N; //Number of neurons
    int i, j, k, counterr, counterV; // counters
    double aux, h, /*time increment*/ t, /*time*/tmax, tmin;
    double Vp, /*peak*/ tau, average, averageV, averager, /*average potential*/ rate; /*firing rate*/
    double Vthres, /*Vthreshold*/ Vr; /*Vrest*/
    int sum; // Number of neurons for the average of the potential.
    double tpot[10000]; // time when the potential was sent
    double eta[10000];
    double s;//mean synaptic activation
    double etamedia, dt, /*time we use for rate*/ Ndt; //N*dt
    double I[10000]; //Input current
    double V[10000]; //Neuron potential
    double J; //Synaptic weight
    double Js; //J*s
    const double pi = atan(1)*4; //PI
    double rast[300]; //Raster vector of 300 neurons;
    bool spike[10000]; //spiking neurons
    double refrac_period[10000]; //Refractory period for every neuron.
    double taverage;
    int conectivity [10000]; //Number of conexions of every Neuron;
    vector< vector<int> > conexion; //Vector de conexiones

    cout<<"eta"<<"  "<<"   V  "<<"  "<<"   r"<<endl;
    N = 10000;
    J = 15.0;
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
    ofstream file1;
    ofstream file2;

     //Generador de aleatorios
    random_device rd;
    mt19937 gen(45454536345331);
    uniform_real_distribution<> dis(0, 100);

    file1.open("archivostxt/meanaverageV_J15_ro1_Vo0-2.txt");
    file2.open("archivostxt/meanrate_J15_ro1_Vo0-2.txt");


    //Creamos el número de filas del vector
    for (i = 0; i<N; i++)
    {
        conexion.push_back(vector<int>());
    }

    //Creamos las distintas conectividades
    for (i=0;i<N;i++)
        {
            for(j=i+1;j<N;j++)
            {
               if (dis(gen)<conectividad)
               {
                    conexion[i].push_back(j);
                    conexion[j].push_back(i);
                    cout<<i<<"  "<<j<<"    ";
               }

            }
            conectivity[i]= conexion[i].size()
        }

    //Inicializamos el periodo refractario. Luego cambiaremos este valor
    for (j=0;j<N;j++)
    {
        refrac_period[j]= 2.0/Vp;
    }


    for (etamedia = -6; etamedia <=-4; etamedia+=0.2)
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
            V[j]= 0;
        }

        // No Neuron has been fired yet
        for (j=0;j<N; j++)
        {
            tpot[j]=-10000;
        }

        rate = 0;

            // We calculate the potential of every Neuron using Gauss method for differential equations. Also
            // we calculate the average potential.
            average = 0;
            sum = 0;

            // We calculate the rate each 10^-2 s
            int ftime = t*10000;

            for (j=0;j <N; j++)
            {
                // We calculate s for Neurona j
                for (i = 0; i < conectivity[j]; i++)
                {
                    neurona_conectada = conexion[j][i];
                    if((t - tpot[neurona_conectada] <= tau)&&(t - tpot[neurona_conectada] >= 0))
                    {
                        s += 1;
                    }
                }

                s = 1.0*s/(conectivity[j]*1.0*tau);
                Js = J*s;

                spike[j] = false;
                //Only if the neuron is not in refractory period.
                if (t - tpot[j] > refrac_period[j])
                {
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
                if (t > taverage)
                {
                    rate = 1.0*rate/(Ndt);
                    averager +=rate;
                    counterr += 1;
                }
                rate = 0;
            }

        }

        file1<<etamedia<<"  "<<averageV/counterV<<endl;
        file2<<etamedia<<"  "<<averager/counterr<<endl;
        cout<<etamedia<<"  "<<averageV/counterV<<"  "<<averager/counterr<<endl;
    }
    file1.close();
    file2.close();
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



