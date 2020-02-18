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
double current(double t);
void Euler(double &V, double I, double h);

int main(void)
{
    int N; //Number of neurons
    int i, j, k; // counters
    int sum; // Number of neurons for the average of the potential.
    int neuron_connected; //number of the connected neuron
    int connectivity [10000]; //Number of conexions of every Neuron;
    double nconex;
    double J; //Synaptic weight
    double Js; //J*s
    double aux, h, /*time increment*/ t, /*time*/tmax, tmin, s;//mean synaptic activation
    double Vp, /*peak*/ tau, average, /*average potential*/ rate; /*firing rate*/
    double Vthres, /*Vthreshold*/ Vr; /*Vrest*/
    double etamedia, dt, /*time we use for rate*/ Ndt; //N*dt
    double prob_connection; //probability for two neurons to be connected
    double tpot[10000]; // time when the potential was sent
    double eta[10000];
    double I[10000]; //Input current
    double V[10000]; //Neuron potential
    double rast[300]; //Raster vector of 300 neurons;
    double refrac_period[10000]; //Refractory period for every neuron.
    const double pi = atan(1)*4; //PI
    bool spike[10001]; //spiking neurons
    vector< vector<int> > connection; //Vector of connections

    //Generador de aleatorios
    random_device rd;
    mt19937 gen(45454536345331);
    uniform_real_distribution<> dis(-102, 98);
    uniform_real_distribution<> con(0,100);
    uniform_int_distribution <> disint(1,10000);

    N = 10000;
    J = 15.0;
    tau = 1.0e-3;
    h = 1.0e-4;
    tmax = 40;
    tmin = -10;
    Vp = 100.0;
    Vthres = Vp;
    Vr = -Vp;
    etamedia = -5.0;
    rate = 0;
    dt = 1e-2;
    Ndt = N*dt;
    prob_connection = 1;

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
                    nconex += 1;
               }

            }
            connectivity[i]= connection[i].size();
        }
        cout << 1.0*nconex/N;

    // We generate a provisional refractory period for every Neuron
    for (j = 0 ;j < N;j++)
    {
        refrac_period[j]= 2.0/Vp;
    }

    // We calculate eta for each Neuron. Only needed once.
    for (j =0; j < N; j++)
    {
        aux = (pi/2.0)*((2.0*(j+1)-N-1)/(N+1.0));
        eta[j] = etamedia + tan(aux);
    }

    // V are distributed randomly
    for (j = 0;j < N; j++)
    {
        V[j]= Vr;
    }

    // No Neuron has been fired yet
    for (j=0;j < N; j++)
    {
        tpot[j]=-10000;
    }



    ofstream file1;
    ofstream file2;
    ofstream file3;
    ofstream file5;
    ofstream file6;


    file1.open("archivostxt/potunaneurona.txt");
    file2.open("archivostxt/meanpotential.txt");
    file3.open("archivostxt/firingrate.txt");
    file5.open("archivostxt/raster.txt");
    file6.open("archivostxt/s.txt");


    rate = 0;

    // We select Neurons for the raster
    for (i=0; i<300; i++)
    {
        rast[i] = disint(gen);
    }



    for (t=tmin; t <= tmax; t += h)

    {
        // We calculate the potential of every Neuron using Gauss method for differential equations. Also
        // we calculate the average potential.
        average = 0;
        sum = 0;

        // We calculate the rate each 10^-2 s
        int ftime = t*10000;

        for (j=0;j < N; j++)
        {
            s = 0;
            spike[j] = false;
            //Only if the neuron is not in refractory period.
            if (t - tpot[j] > refrac_period[j])
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
                Js = J*s;

                I[j] = eta[j] + Js + current(t);
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


        file1 << t <<"  "<< V[5000] << "  " << I[5000] <<endl;

        file2 << t <<"  "<< average/(1.0*sum) << endl;
        if (ftime%100 == 0)
        {
            file3 << t <<"  "<< 1.0*rate/(Ndt) << endl;
            rate = 0;
        }

        // Finally we plot the raster
        for ( i=0;i<300;i++)
        {
            k = rast[i];
            if (spike[k] == true)
            {
                file5 <<t<<"  "<<i<<endl;
            }
        }

    }

    file1.close();
    file2.close();
    file3.close();
    file5.close();
    file6.close();

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

double current(double t)
//External current
{
     double I;
     if ((t<0)||(t>30))
     {
         I = 0.0;
     }
     else
     {
         I = 3.0;
     }
     return I;
//    double I;
//    const double pi = atan(1)*4; //PI
//    I = 3*sin(pi/20*t);
//    return I;
}

void Euler(double &V, double I, double h)
    {
        V = V + h*dVdt(V,I);
        return;
    }

