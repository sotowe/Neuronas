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


int main(void)
{
    int N = 10;
    int i, j, tamano;
    int conectividad;
    vector< vector<int> > conexion(N);

    conectividad = 50;

     //Generador de aleatorios
    random_device rd;
    mt19937 gen(45454536345331);
    uniform_real_distribution<> dis(0, 100);

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
            cout <<endl;
        }
        cout<<endl;
        for (i = 0; i<N; i++)
        {
           tamano = conexion[i].size();
           for(j = 0; j<tamano; j++)
           {
               cout<<conexion[i][j]<<"  ";
           }
           cout<<endl;
        }


    return(0);
}
