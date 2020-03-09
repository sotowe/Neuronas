//Libraries
#include<iostream>
#include<cstdlib>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include<vector>
#include <algorithm>

using namespace std;


// FUNCTIONS
double dVdt(double V, double I);
double current(double t);
void Euler(double &V, double I, double h);

int main(void)
{
    int N; //Number of neurons
    int i, j, k, l, m; // counters
    int neuron_connected; //number of the connected neuron
    int Nraster; //Number of neurons for the raster.
    int connectivity [10000]; //Number of conexions of every Neuron;
    int sum[5]; // Number of neurons for the average of the potential.
    double nconex;
    double J; //Synaptic weight
    double Js; //J*s
    double aux, h, /*time increment*/ t, /*time*/tmax, tmin, s;//mean synaptic activation
    double Vp, /*peak*/ tau;
    double Vthres, /*Vthreshold*/ Vr; /*Vrest*/
    double etamedia, dt, /*time we use for rate*/ Ndt, /*N*dt*/N_4dt; //N_4*dt
    double prob_connection1; //probability for two neurons to be connected in one module
    double prob_connection2; //probability for two neurons to be connected between different modules
    double tpot[10000]; // time when the potential was sent
    double average [5], /*average potential*/ rate[5]; /*firing rate*/
    double eta[10000];
    double I[10000]; //Input current
    double V[10000]; //Neuron potential
    double rast[400]; //Raster vector of 400 neurons;
    double refrac_period[10000]; //Refractory period for every neuron.
    const double pi = atan(1)*4; //PI
    bool spike[10000]; //spiking neurons
    vector< vector<int> > connection; //Vector of connections

    //Generador de aleatorios
    random_device rd;
    mt19937 gen(45454536345331);
    uniform_real_distribution<> dis(-102, 98);
    uniform_real_distribution<> con(0,100);
    uniform_int_distribution <> disint(0,2499);

    N = 10000;
    Nraster = 400;
    J = 15.0;
    tau = 1.0e-3;
    h = 1.0e-4;
    tmax = 80;
    tmin = -10;
    Vp = 100.0;
    Vthres = Vp;
    Vr = -Vp;
    etamedia = -5.0;
    for (k=0;k<5;k++) rate[k] = 0;
    dt = 1e-2;
    Ndt = N*dt;
    prob_connection1 = 0.4;
    prob_connection2 = 0.04;

    ofstream file1;
    ofstream file2;
    ofstream file3;
    ofstream file4;
    file1.open("archivostxt/pot_modular.txt");
    file2.open("archivostxt/rate_modular.txt");
    file3.open("archivostxt/raster_modular.txt");
    file4.open("archivostxt/current.txt");

    //Creamos el número de filas del vector
    for (i = 0; i<N; i++)
    {
        connection.push_back(vector<int>());
    }
    int N_4 = N/4;
    N_4dt = N_4*dt;
    //Creamos las distintas conectividades intramodulares
    for (k=0; k<4; k++)
    {
        nconex = 0;
        for (i=N_4*k;i<N_4*(k+1);i++)
        {
            for(j=i+1;j<N_4*(k+1);j++)
            {
               if (con(gen)<prob_connection1)
               {
                    connection[i].push_back(j);
                    connection[j].push_back(i);
                    nconex += 1;
               }
            }
        }
        cout << 2.0*nconex/N_4<<endl;
    }
    //Creamos las distintas conectividades intermodulares
    for (k=0; k<3; k++)
    {
        nconex = 0;
        for (i=N_4*k;i<N_4*(k+1);i++)
        {
            for(j=i+N_4;j<N_4*(k+2);j++)
            {
               if (con(gen)<prob_connection2)
               {
                    connection[i].push_back(j);
                    connection[j].push_back(i);
                    nconex += 1;
               }
            }
        }
        cout << 2.0*nconex/N_4<<endl;
    }
    nconex = 0;
    for (i=0;i<N_4;i++)
    {
        for(j=i+3*N_4;j<N;j++)
        {
           if (con(gen)<prob_connection2)
           {
                connection[i].push_back(j);
                connection[j].push_back(i);
                nconex += 1;
           }
        }
    }
   cout << 2.0*nconex/N_4<<endl;
   for (i = 0; i<N; i++)
   {
       connectivity[i]= connection[i].size();
   }

    // Generamos un periodo refractario provisional para todas las neuronas.
    for (j = 0 ;j < N;j++)
    {
        refrac_period[j]= 1.0/Vp;
    }

    // We calculate eta for each Neuron. We distribute equally eta between the different modules.
    for (j =0; j < N_4; j++)
    {
        aux = (pi/2.0)*((2.0*(j+1)-N_4-1)/(N_4+1.0));
        eta[j] = etamedia + tan(aux);
        for (k=1; k<4; k++)
        {
            eta[j + k*N_4]=eta[j];
        }
    }

    // V todos en reposo
    for (j = 0;j < N; j++)
    {
        V[j]= Vr;
    }

    // No Neuron has been fired yet
    for (j=0;j < N; j++)
    {
        tpot[j]=-10000;
    }

    // Seleccionamos las neuronas para el raster. Estarán en orden según módulo y ordenadas dentro de él.
    for (k = 0; k<4; k++)
    {
        for (i=100*k; i<100*(k+1); i++)
        {
            rast[i] = disint(gen)+2500*k;
        }
    }

    for (t=tmin; t <= tmax; t += h)
    {
        file4<<t<<"  "<<current(t)<<endl;
        // We calculate the potential of every Neuron using Gauss method for differential equations. Also
        // we calculate the average potential. Y el potencial medio de cada módulo
        for (i=0; i<5; i++)
        {
            average[i]=0;
            sum[i]=0;
        }
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

                // Solo vamos a enviar corriente a uno de los módulos
                if (j<N_4)
                {
                    I[j] = eta[j] + Js + current(t);
                }
                else
                {
                    I[j] = eta[j] + Js;
                }

                Euler(V[j],I[j],h);

                 // If the potential of the neuron reach Vthres, the neuron sends an impulse and goes to a refractory state.
                if (V[j]>= Vthres)
                {
                    refrac_period[j]=1.0/V[j];
                    tpot[j] = t + 1.0/V[j];
                    V[j] = -V[j];
                    rate[0] ++;
                    if (j < N_4) rate[1]++;
                    else if (j < 2*N_4) rate[2]++;
                    else if (j < 3*N_4) rate[3]++;
                    else rate[4]++;
                    spike[j] = true;
                }
            }
            // We calculate the average potential in the network. Only no refractory.
            if ((t - tpot[j] > refrac_period[j])&&(V[j]==V[j]))
            {
                average[0] += V[j];
                sum[0] ++;

                if (j < N_4)
                {
                    average[1] += V[j];
                    sum[1] ++;
                }
                else if (j < 2*N_4)
                {
                    average[2] += V[j];
                    sum[2] ++;
                }
                else if (j < 3*N_4)
                {
                    average[3] += V[j];
                    sum[3] ++;
                }
                else
                {
                    average[4] += V[j];
                    sum[4] ++;
                }
            }
        }
        file1 << t <<"  ";
        for (k=0; k<=4; k++)
        {
            average[k]= average[k]/(1.0*sum[k]);
            file1 << average[k] << "   ";
        }
        file1<<endl;
        if (ftime%100 == 0)
        {
            file2 << t <<"  "<< 1.0*rate[0]/Ndt <<"  ";
            rate[0] = 0;
            for (k=1;k<=4;k++)
            {
                rate[k] = 1.0*rate[k]/N_4dt;
                file2 << rate[k]<< "  ";
                rate[k]= 0;
            }
            file2<<endl;
        }

        // Finally we plot the raster
        for ( i=0;i<Nraster;i++)
        {
            k = rast[i];
            if (spike[k] == true)
            {
                file3 <<t<<"  "<<i<<endl;
            }
        }
    }

    file1.close();
    file2.close();
    file3.close();
    file4.close();

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
//     double I;
//     if ((t<0)||(t>30))
//     {
//         I = 0.0;
//     }
//     else
//     {
//         I = 3.0;
//     }
//     return I;
    double I;
    const double pi = atan(1)*4; //PI
    if (t<0)
    {
        I = 0;
    }
    else
    {
        I = 3*sin(pi/20*t);
    }
    return I;
}

void Euler(double &V, double I, double h)
    {
        V = V + h*dVdt(V,I);
        return;
    }

