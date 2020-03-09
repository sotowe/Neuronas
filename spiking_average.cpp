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
    double eta_ini, eta_fin;
    double tpot[10000]; // time when the potential was sent
    double eta[10000];
    double I[10000]; //Input current
    double V[10000]; //Neuron potential
    double refrac_period[10000]; //Refractory period for every neuron.
    const double pi = atan(1)*4; //PI
    bool spike[10000]; //spiking neurons


    cout<<"eta"<<"  "<<"   V  "<<"  "<<"   r"<<"    "<<"errorV"<<"   "<<"errorr"<<endl;
    N = 10000;
    J = 20.0;
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
    ofstream file1;
    ofstream file2;
    ofstream file3;

    for (J = 15; J <21; J=J+5)
    {
        if (J == 5)
        {
            file1.open("archivostxt/meanaverageV_J05_ro0_Vo-3.txt");
            file2.open("archivostxt/meanrate_J05_ro0_Vo-3.txt");
            cout<<J<<endl;
        }
        else if (J==10)
        {
            file1.open("archivostxt/meanaverageV_J10_ro0_Vo-3.txt");
            file2.open("archivostxt/meanrate_J10_ro0_Vo-3.txt");
            eta_ini = -3;
            eta_fin = -2;
             cout<<J<<endl;
        }
        else if (J==15)
        {
            file1.open("archivostxt/meanaverageV_J15_ro0_Vo-3.txt");
            file2.open("archivostxt/meanrate_J15_ro0_Vo-3.txt");
            eta_ini = -4;
            eta_fin = -3;
             cout<<J<<endl;
        }
        else if (J==20)
        {
            file1.open("archivostxt/meanaverageV_J20_ro0_Vo-3.txt");
            file2.open("archivostxt/meanrate_J20_ro0_Vo-3.txt");
            eta_ini = -4;
            eta_fin = -3;
             cout<<J<<endl;
        }


        for (etamedia = -11; etamedia <= 0; etamedia++)
        {
            for (j=0;j < N;j++)
            {
                refrac_period[j]= 1.0/Vp;
            }
            if (J == 5)
            {
                if (etamedia == -11) file3.open("archivostxt/rate_11_J05_r0.txt");
               else if (etamedia == -10) file3.open("archivostxt/rate_J05_10_r0.txt");
               else if (etamedia == -9) file3.open("archivostxt/rate_9_J05_r0.txt");
               else if (etamedia == -8) file3.open("archivostxt/rate_8_J05_r0.txt");
               else if (etamedia == -7) file3.open("archivostxt/rate_7_J05_r0.txt");
               else if (etamedia == -6) file3.open("archivostxt/rate_6_J05_r0.txt");
               else if (etamedia == -5) file3.open("archivostxt/rate_5_J05_r0.txt");
               else if (etamedia == -4) file3.open("archivostxt/rate_4_J05_r0.txt");
               else if (etamedia == -3) file3.open("archivostxt/rate_3_J05_r0.txt");
               else if (etamedia == -2) file3.open("archivostxt/rate_2_J05_r0.txt");
               else if (etamedia == -1) file3.open("archivostxt/rate_1_J05_r0.txt");
               else if (etamedia == 0) file3.open("archivostxt/rate_0_J05_r0.txt");
            }
            else if (J == 10)
            {
                if (etamedia == -11) file3.open("archivostxt/rate_11_J10_r0.txt");
               else if (etamedia == -10) file3.open("archivostxt/rate_J10_10_r0.txt");
               else if (etamedia == -9) file3.open("archivostxt/rate_9_J10_r0.txt");
               else if (etamedia == -8) file3.open("archivostxt/rate_8_J10_r0.txt");
               else if (etamedia == -7) file3.open("archivostxt/rate_7_J10_r0.txt");
               else if (etamedia == -6) file3.open("archivostxt/rate_6_J10_r0.txt");
               else if (etamedia == -5) file3.open("archivostxt/rate_5_J10_r0.txt");
               else if (etamedia == -4) file3.open("archivostxt/rate_4_J10_r0.txt");
               else if (etamedia == -3) file3.open("archivostxt/rate_3_J10_r0.txt");
               else if (etamedia == -2) file3.open("archivostxt/rate_2_J10_r0.txt");
               else if (etamedia == -1) file3.open("archivostxt/rate_1_J10_r0.txt");
               else if (etamedia == 0) file3.open("archivostxt/rate_0_J10_r0.txt");
            }
             else if (J == 15)
            {
                if (etamedia == -11) file3.open("archivostxt/rate_11_J15_r0.txt");
               else if (etamedia == -10) file3.open("archivostxt/rate_J15_10_r0.txt");
               else if (etamedia == -9) file3.open("archivostxt/rate_9_J15_r0.txt");
               else if (etamedia == -8) file3.open("archivostxt/rate_8_J15_r0.txt");
               else if (etamedia == -7) file3.open("archivostxt/rate_7_J15_r0.txt");
               else if (etamedia == -6) file3.open("archivostxt/rate_6_J15_r0.txt");
               else if (etamedia == -5) file3.open("archivostxt/rate_5_J15_r0.txt");
               else if (etamedia == -4) file3.open("archivostxt/rate_4_J15_r0.txt");
               else if (etamedia == -3) file3.open("archivostxt/rate_3_J15_r0.txt");
               else if (etamedia == -2) file3.open("archivostxt/rate_2_J15_r0.txt");
               else if (etamedia == -1) file3.open("archivostxt/rate_1_J15_r0.txt");
               else if (etamedia == 0) file3.open("archivostxt/rate_0_J15_r0.txt");
            }
             else if (J == 20)
            {
                if (etamedia == -11) file3.open("archivostxt/rate_11_J20_r0.txt");
               else if (etamedia == -10) file3.open("archivostxt/rate_J20_10_r0.txt");
               else if (etamedia == -9) file3.open("archivostxt/rate_9_J20_r0.txt");
               else if (etamedia == -8) file3.open("archivostxt/rate_8_J20_r0.txt");
               else if (etamedia == -7) file3.open("archivostxt/rate_7_J20_r0.txt");
               else if (etamedia == -6) file3.open("archivostxt/rate_6_J20_r0.txt");
               else if (etamedia == -5) file3.open("archivostxt/rate_5_J20_r0.txt");
               else if (etamedia == -4) file3.open("archivostxt/rate_4_J20_r0.txt");
               else if (etamedia == -3) file3.open("archivostxt/rate_3_J20_r0.txt");
               else if (etamedia == -2) file3.open("archivostxt/rate_2_J20_r0.txt");
               else if (etamedia == -1) file3.open("archivostxt/rate_1_J20_r0.txt");
               else if (etamedia == 0) file3.open("archivostxt/rate_0_J20_r0.txt");
            }
            counterr = 0;
            counterV = 0;
            averageV = 0;
            averager = 0;
            averageV2 = 0;
            averager2 = 0;

            // We calculate eta for each Neuron. Only needed once.
            for (j = 0; j <  N; j++)
            {
                aux = (pi/2.0)*((2.0*(j+1)-N-1)/(N+1.0));
                eta[j] = etamedia + tan(aux);
            }

            // V are all 0
            for (j=0;j<N; j++)
            {
                V[j]= -3;
            }

            // No Neuron has been fired yet
            for (j=0;j<N; j++)
            {
                tpot[j]=-10000;
            }

            rate = 0;
            for (t=tmin; t <= tmax; t += h)
            {
             // We calculate  the mean synaptic activation for every Neuron. The first 15 seconds we have r = 0;
                if (t<0)
                {
                    s = 0;
                }
                else
                {
                    s = 0;
                    for (j=0; j<N; j++)
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

                for (j=0;j < N; j++)
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
                            refrac_period[j]=1.0/V[j];
                            tpot[j] = t + 1.0/V[j];
                            V[j] = -V[j];
                            rate = rate += 1;
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
            file3.close();
        }
        file1.close();
        file2.close();

        if (J == 5)
        {
            file1.open("archivostxt/meanaverageV_J05_ro1_Vo0.txt");
            file2.open("archivostxt/meanrate_J05_ro1_Vo0.txt");
            cout<<J<<endl;
        }
        else if (J==10)
        {
            file1.open("archivostxt/meanaverageV_J10_ro1_Vo0.txt");
            file2.open("archivostxt/meanrate_J10_ro1_Vo0.txt");
            eta_ini = -3;
            eta_fin = -2;
             cout<<J<<endl;
        }
        else if (J==15)
        {
            file1.open("archivostxt/meanaverageV_J15_ro1_Vo0.txt");
            file2.open("archivostxt/meanrate_J15_ro1_Vo0.txt");
            eta_ini = -6;
            eta_fin = -5;
             cout<<J<<endl;
        }
        else if (J==20)
        {
            file1.open("archivostxt/meanaverageV_J20_ro1_Vo0.txt");
            file2.open("archivostxt/meanrate_J20_ro1_Vo0.txt");
            eta_ini = -11;
            eta_fin = -10;
             cout<<J<<endl;
        }

        for (j=0;j< N;j++)
        {
            refrac_period[j]= 1.0/Vp;
        }
        for (etamedia = -11; etamedia <=0; etamedia++)
        {
            if (J == 5)
            {
                if (etamedia == -11) file3.open("archivostxt/rate_11_J05_r1.txt");
               else if (etamedia == -10) file3.open("archivostxt/rate_J05_10_r1.txt");
               else if (etamedia == -9) file3.open("archivostxt/rate_9_J05_r1.txt");
               else if (etamedia == -8) file3.open("archivostxt/rate_8_J05_r1.txt");
               else if (etamedia == -7) file3.open("archivostxt/rate_7_J05_r1.txt");
               else if (etamedia == -6) file3.open("archivostxt/rate_6_J05_r1.txt");
               else if (etamedia == -5) file3.open("archivostxt/rate_5_J05_r1.txt");
               else if (etamedia == -4) file3.open("archivostxt/rate_4_J05_r1.txt");
               else if (etamedia == -3) file3.open("archivostxt/rate_3_J05_r1.txt");
               else if (etamedia == -2) file3.open("archivostxt/rate_2_J05_r1.txt");
               else if (etamedia == -1) file3.open("archivostxt/rate_1_J05_r1.txt");
               else if (etamedia == 0) file3.open("archivostxt/rate_0_J05_r1.txt");
            }
            else if (J == 10)
            {
                if (etamedia == -11) file3.open("archivostxt/rate_11_J10_r1.txt");
               else if (etamedia == -10) file3.open("archivostxt/rate_J10_10_r1.txt");
               else if (etamedia == -9) file3.open("archivostxt/rate_9_J10_r1.txt");
               else if (etamedia == -8) file3.open("archivostxt/rate_8_J10_r1.txt");
               else if (etamedia == -7) file3.open("archivostxt/rate_7_J10_r1.txt");
               else if (etamedia == -6) file3.open("archivostxt/rate_6_J10_r1.txt");
               else if (etamedia == -5) file3.open("archivostxt/rate_5_J10_r1.txt");
               else if (etamedia == -4) file3.open("archivostxt/rate_4_J10_r1.txt");
               else if (etamedia == -3) file3.open("archivostxt/rate_3_J10_r1.txt");
               else if (etamedia == -2) file3.open("archivostxt/rate_2_J10_r1.txt");
               else if (etamedia == -1) file3.open("archivostxt/rate_1_J10_r1.txt");
               else if (etamedia == 0) file3.open("archivostxt/rate_0_J10_r1.txt");
            }

             else if (J == 15)
            {
                if (etamedia == -11) file3.open("archivostxt/rate_11_J15_r1.txt");
               else if (etamedia == -10) file3.open("archivostxt/rate_J15_10_r1.txt");
               else if (etamedia == -9) file3.open("archivostxt/rate_9_J15_r1.txt");
               else if (etamedia == -8) file3.open("archivostxt/rate_8_J15_r1.txt");
               else if (etamedia == -7) file3.open("archivostxt/rate_7_J15_r1.txt");
               else if (etamedia == -6) file3.open("archivostxt/rate_6_J15_r1.txt");
               else if (etamedia == -5) file3.open("archivostxt/rate_5_J15_r1.txt");
               else if (etamedia == -4) file3.open("archivostxt/rate_4_J15_r1.txt");
               else if (etamedia == -3) file3.open("archivostxt/rate_3_J15_r1.txt");
               else if (etamedia == -2) file3.open("archivostxt/rate_2_J15_r1.txt");
               else if (etamedia == -1) file3.open("archivostxt/rate_1_J15_r1.txt");
               else if (etamedia == 0) file3.open("archivostxt/rate_0_J15_r1.txt");
            }
             else if (J == 20)
            {
                if (etamedia == -11) file3.open("archivostxt/rate_11_J20_r1.txt");
               else if (etamedia == -10) file3.open("archivostxt/rate_J20_10_r1.txt");
               else if (etamedia == -9) file3.open("archivostxt/rate_9_J20_r1.txt");
               else if (etamedia == -8) file3.open("archivostxt/rate_8_J20_r1.txt");
               else if (etamedia == -7) file3.open("archivostxt/rate_7_J20_r1.txt");
               else if (etamedia == -6) file3.open("archivostxt/rate_6_J20_r1.txt");
               else if (etamedia == -5) file3.open("archivostxt/rate_5_J20_r1.txt");
               else if (etamedia == -4) file3.open("archivostxt/rate_4_J20_r1.txt");
               else if (etamedia == -3) file3.open("archivostxt/rate_3_J20_r1.txt");
               else if (etamedia == -2) file3.open("archivostxt/rate_2_J20_r1.txt");
               else if (etamedia == -1) file3.open("archivostxt/rate_1_J20_r1.txt");
               else if (etamedia == 0) file3.open("archivostxt/rate_0_J20_r1.txt");
            }

            counterr = 0;
            counterV = 0;
            averageV = 0;
            averager = 0;
            averageV2 = 0;
            averager2 = 0;

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
            for (t=tmin; t <= tmax; t += h)
            {
             // We calculate  the mean synaptic activation for every Neuron. The first 5 seconds we have r = 1;
                    if (t<0)
                    {
                        s = 1;
                    }
                    else
                    {
                        s = 0;
                        for (j=0; j<N; j++)
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

                for (j=0;j <N; j++)
                {
                    //Only if the neuron is not in refractory period.
                    if (t - tpot[j] > refrac_period[j])
                    {

                        I[j] = eta[j] + Js;
                        Euler(V[j],I[j],h);

                         // If the potential of the neuron reach Vthres, the neuron sends an impulse and
                         //goes to a refractory state.
                        if (V[j]>= Vthres)
                        {
                            refrac_period[j]=1.0/V[j];
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
            file3.close();
        }
        file1.close();
        file2.close();
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
    error = 1.96*sqrt(1.0*(media2 - medidas*media*media))/medidas;
    return error;
}
