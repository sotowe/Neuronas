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
    int connectivity [10000]; //Number of conexions of every Neuron;
    int neuron_connected; //number of the connected neuron
    int conditions; //Preliminary conditions
    double aux, h, /*time increment*/ t, /*time*/tmax, tmin, s;//mean synaptic activation
    double Vp, /*peak*/ tau, average, averageV, averager, /*average potential*/ rate; /*firing rate*/
    double Vthres, /*Vthreshold*/ Vr; /*Vrest*/
    double etamedia, dt, /*time we use for rate*/ Ndt; //N*dt
    double J; //Synaptic weight
    double Js; //J*s
    double taverage;
    double prob_connection; //probability for two neurons to be connected
    double errorr, errorV;
    double averageV2, averager2; //sum(V_i^2) and sum(r_i^2)
    double nconex; //Number of connections
    double prob_connection1; //probability for two neurons to be connected in one module
    double prob_connection2; //probability for two neurons to be connected between different modules
    double eta_ini, eta_fin;
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

    cout<<"eta"<<"  "<<"   V  "<<"  "<<"   r"<<"    "<<"errorV"<<"   "<<"errorr"<<endl;
    N = 10000;
    J = 5.0;
    tau = 1.0e-3;
    h = 1.0e-4;
    tmax = 40;
    tmin = -15;
    taverage = 20;
    Vp = 100.0;
    Vthres = Vp;
    Vr = -Vp;
    rate = 0;
    dt = 1e-2;
    Ndt = N*dt;
    prob_connection1 = 0.4;
    prob_connection2 = 0.04;

    ofstream file1;
    ofstream file2;
    ofstream file3;

    //Creamos el número de filas del vector
    for (i = 0; i<N; i++)
    {
        connection.push_back(vector<int>());
    }
    int N_4 = N/4;
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

    // We generate a provisional refractory period for every Neuron

    for (J = 10; J <= 21; J=J+5)
    {
        for (conditions = 1; conditions <=2; conditions ++)
        {
            cout<<J<<endl;
            if (J == 5)
            {
                if (conditions == 1)
                {
                    file1.open("archivostxt/meanaverageV_J5_ro0_Vo-3_mod10_01_2.txt");
                    file2.open("archivostxt/meanrate_J5_ro0_Vo-3_mod10_01_2.txt");
                }
                else if (conditions == 2)
                {
                    file1.open("archivostxt/meanaverageV_J5_ro1_Vo0_mod10_01_2.txt");
                    file2.open("archivostxt/meanrate_J5_ro1_Vo0_mod10_01_2.txt");
                }
            }
            else if (J == 10)
            {
                if (conditions == 1)
                {
                    file1.open("archivostxt/meanaverageV_J10_ro0_Vo-3_mod10_01_2.txt");
                    file2.open("archivostxt/meanrate_J10_ro0_Vo-3_mod10_01_2.txt");
                    eta_ini = -3;
                    eta_fin = -2;
                }
                else if (conditions == 2)
                {
                    file1.open("archivostxt/meanaverageV_J10_ro1_Vo0_mod10_01_2.txt");
                    file2.open("archivostxt/meanrate_J10_ro1_Vo0_mod10_01_2.txt");
                    eta_ini = -3;
                    eta_fin = -2;
                }
            }
            else if (J == 15)
            {
                 if (conditions == 1)
                {
                    file1.open("archivostxt/meanaverageV_J15_ro0_Vo-3_mod10_01_2.txt");
                    file2.open("archivostxt/meanrate_J15_ro0_Vo-3_mod10_01_2.txt");
                    eta_ini = -5;
                    eta_fin = -4;
                }
                else if (conditions == 2)
                {
                    file1.open("archivostxt/meanaverageV_J15_ro1_Vo0_mod10_01_2.txt");
                    file2.open("archivostxt/meanrate_J15_ro1_Vo0_mod10_01_2.txt");
                    eta_ini = -7;
                    eta_fin = -6;
                }
            }
            else if (J == 20)
            {
                 if (conditions == 1)
                {
                    file1.open("archivostxt/meanaverageV_J20_ro0_Vo-3_mod10_01_2.txt");
                    file2.open("archivostxt/meanrate_J20_ro0_Vo-3_mod10_01_2.txt");
                    eta_ini = -7;
                    eta_fin = -6;
                }
                else if (conditions == 2)
                {
                    file1.open("archivostxt/meanaverageV_J20_ro1_Vo0_mod10_01_2.txt");
                    file2.open("archivostxt/meanrate_J20_ro1_Vo0_mod10_01_2.txt");
                    eta_ini = -12;
                    eta_fin = -11;
                }
            }
            for (etamedia = eta_ini; etamedia <= eta_fin; etamedia+=0.2)
            {
                for (j=0;j<N;j++)
                {
                    refrac_period[j]= 1.0/Vp;
                }

                counterr = 0;
                counterV = 0;
                averageV = 0;
                averager = 0;
                averageV2 = 0;
                averager2 = 0;

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

                if (conditions ==1)
                {
                    // V are all 0
                    for (j=0;j<N; j++)
                    {
                        V[j]= -3.0;
                    }
                }
                else if (conditions ==2)
                {
                    // V are all 0
                    for (j=0;j<N; j++)
                    {
                        V[j]= 0.0;
                    }
                }

                // No Neuron has been fired yet
                for (j=0;j<N; j++)
                {
                    tpot[j]=-10000;
                }

                rate = 0;
                for (t=tmin; t <= tmax; t += h)
                {

                    // We calculate the rate each 10^-2 s
                    int ftime = t*10000;
                    average = 0;
                    sum = 0;
                    for (j=0;j<N; j++)
                    {
                        s = 0;
                        //Only if the neuron is not in refractory period.
                        if (t - tpot[j] > refrac_period[j])
                        {
                            if (t<0)
                            {
                                if (conditions==1)
                                {
                                    s = 0;
                                }
                                else if (conditions == 2)
                                {
                                    s = 1;
                                }

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
                                refrac_period[j]=1.0/V[j];
                                tpot[j] = t + 1.0/V[j];
                                V[j] = -V[j];
                                rate = rate += 1;
                            }

                        }
                        // We calculate the average potential in the network. Only no refractory.
                        if ((t - tpot[j] > refrac_period[j])&&(V[j]==V[j]))
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
                file3.close();
            }
            file1.close();
            file2.close();

        }
    }

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
    error = sqrt(1.0*(media2 - medidas*media*media))/medidas;
    return error;
}