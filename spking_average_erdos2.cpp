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
    double eta_ini, eta_fin;
    double tpot[10000]; // time when the potential was sent
    double eta[10000];
    double refrac_period[10000]; //Refractory period for every neuron.
    double I[10000]; //Input current
    double V[10000]; //Neuron potential
    const double pi = atan(1)*4; //PI
    vector< vector<int> > connection; //Vector of connections

    //Generador de aleatorios
    random_device rd;
    mt19937 gen(45454536345331);
    uniform_real_distribution<> con(0,100);

    cout<<"eta"<<"  "<<"   V  "<<"  "<<"   r"<<"    "<<"errorV"<<"   "<<"errorr"<<endl;
    N = 10000;
    tau = 1.0e-3;
    h = 1.0e-4;
    tmax = 60;
    tmin = -15;
    taverage = 30;
    Vp = 100.0;
    Vthres = Vp;
    Vr = -Vp;
    rate = 0;
    dt = 1e-2;
    Ndt = N*dt;
    prob_connection = 0.1;

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

    ofstream file1;
    ofstream file2;
    ofstream file3;
    file3.open("archivostxt/rate4,6_J15_ro_con10.txt");
    for (J = 15; J <= 15; J=J+5)
    {
        for (conditions = 1; conditions <2; conditions ++)
        {
            cout<<J<<endl;
            if (J == 5)
            {
                 if (conditions == 1)
                {
                    file1.open("archivostxt/meanaverageV_J5_ro0_Vo-3_con10_2.txt");
                    file2.open("archivostxt/meanrate_J5_ro0_Vo-3_con10_2.txt");
                }
                else if (conditions == 2)
                {
                    file1.open("archivostxt/meanaverageV_J5_ro1_Vo0_con10_2.txt");
                    file2.open("archivostxt/meanrate_J5_ro1_Vo0_con10_2.txt");
                }
            }
            else if (J == 10)
            {
                if (conditions == 1)
                {
                    file1.open("archivostxt/meanaverageV_J10_ro0_Vo-3_con10_2.txt");
                    file2.open("archivostxt/meanrate_J10_ro0_Vo-3_con10_2.txt");
                    eta_ini = -3;
                    eta_fin = -2;
                }
                else if (conditions == 2)
                {
                    file1.open("archivostxt/meanaverageV_J10_ro1_Vo0_con10_2.txt");
                    file2.open("archivostxt/meanrate_J10_ro1_Vo0_con10_2.txt");
                    eta_ini = -3;
                    eta_fin = -2;
                }
            }
            else if (J == 15)
            {
                 if (conditions == 1)
                {
                    file1.open("archivostxt/meanaverageV_J15_ro0_Vo-3_con10_2.txt");
                    file2.open("archivostxt/meanrate_J15_ro0_Vo-3_con10_2.txt");
                    eta_ini = -4.6;
                    eta_fin = -4.5;
                }
                else if (conditions == 2)
                {
                    file1.open("archivostxt/meanaverageV_J15_ro1_Vo0_con10_2.txt");
                    file2.open("archivostxt/meanrate_J15_ro1_Vo0_con10_2.txt");
                    eta_ini = -7;
                    eta_fin = -6;
                }
            }
            else if (J == 20)
            {
                 if (conditions == 1)
                {
                    file1.open("archivostxt/meanaverageV_J20_ro0_Vo-3_con10_2.txt");
                    file2.open("archivostxt/meanrate_J20_ro0_Vo-3_con10_2.txt");
                    eta_ini = -7;
                    eta_fin = -6;
                }
                else if (conditions == 2)
                {
                    file1.open("archivostxt/meanaverageV_J20_ro1_Vo0_con10_2.txt");
                    file2.open("archivostxt/meanrate_J520_ro1_Vo0_con10_2.txt");
                    eta_ini = -12;
                    eta_fin = -11;
                }
            }
            for (etamedia = eta_ini; etamedia <= eta_fin; etamedia+=0.2)
            {
                // We generate a provisional refractory period for every Neuron
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

                // We calculate eta for each Neuron. Just needed once.
                for (j =0; j < N; j++)
                {
                    aux = (pi/2.0)*((2.0*(j+1)-N-1)/(N+1.0));
                    eta[j] = etamedia + tan(aux);
                }

                if (conditions ==1)
                {
                    // V are all -3
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
                    for (j=0;j < N; j++)
                    {
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
                                // Calculamos s para la neurona j
                                s = 0;
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
                        // We calculate the average potential in the network. Only no refractory state.
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
                        file3<<t<<"  "<<rate<<endl;
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
