//Libraries
#include<iostream>
#include<cstdlib>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;


int main(void)
{
    int n;
   random_device rd;
   mt19937 gen(rd());
   uniform_real_distribution<> dis(-1, 1);
   for(n = 1; n<10;n++)
   {
       cout <<dis(gen)<<endl;
   }

}



