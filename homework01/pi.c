#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "common.h"
#include <inttypes.h>
#include <time.h>


#define PI 3.1415926535


void usage(int argc, char** argv);
double calcPi_Serial(int num_steps);
double calcPi_P1(int num_steps);
double calcPi_P2(int num_steps);


int main(int argc, char** argv)
{
    // get input values
    uint32_t num_steps = 100000;
    if(argc > 1) {
        num_steps = atoi(argv[1]);
    } else {
        usage(argc, argv);
        printf("using %"PRIu32"\n", num_steps);
    }
    fprintf(stdout, "The first 10 digits of Pi are %0.10f\n", PI);


    // set up timer
    uint64_t start_t;
    uint64_t end_t;
    InitTSC();

    // generate filename
    int num_cores = atoi(argv[2]);
    char filename[64];
    snprintf(filename, sizeof(filename), "%d-result.txt", num_cores);

    // create file for outputs
    FILE* results;
    results = fopen(filename, "ab+");
    if (results == NULL) {
        printf("Cannot create file");
        return 1;
    }

    // calculate in serial 
    start_t = ReadTSC();
    double Pi0 = calcPi_Serial(num_steps);
    end_t = ReadTSC();
    printf("Time to calculate Pi serially with %"PRIu32" steps is: %g\n",
           num_steps, ElapsedTime(end_t - start_t));
    printf("Pi is %0.10f\n", Pi0);

    fprintf(results, "%g:%0.10f ", ElapsedTime(end_t - start_t), Pi0);
    
    // calculate in parallel with integration
    start_t = ReadTSC();
    double Pi1 = calcPi_P1(num_steps);
    end_t = ReadTSC();
    printf("Time to calculate Pi in // with %"PRIu32" steps is: %g\n",
            num_steps, ElapsedTime(end_t - start_t));
    printf("Pi is %0.10f\n", Pi1);

    fprintf(results, "%g:%0.10f ", ElapsedTime(end_t - start_t), Pi1);



    // calculate in parallel with Monte Carlo
    start_t = ReadTSC();
    double Pi2 = calcPi_P2(num_steps);
    end_t = ReadTSC();
    printf("Time to calculate Pi in // with %"PRIu32" guesses is: %g\n",
           num_steps, ElapsedTime(end_t - start_t));
    printf("Pi is %0.10f\n", Pi2);

    fprintf(results, "%g:%0.10f\n", ElapsedTime(end_t - start_t), Pi2);
    
    return 0;
}


void usage(int argc, char** argv)
{
    fprintf(stdout, "usage: %s <# steps>\n", argv[0]);
}

double calcPi_Serial(int num_steps)
{
    double pi = 0.0;
    int inside_circle = 0;
    double x = 0;
    double y = 0;

    srand(time(NULL));

    for(int i = 0; i < num_steps; i++)
    {
        x = ((double)rand()) / RAND_MAX;
        y = ((double)rand()) / RAND_MAX;
        double coords = (x * x) + (y * y);
        if(coords <= 1)
        {
            inside_circle++;
        }
    }

    pi = 4.0 * ((double)inside_circle / (double)num_steps);
    return pi;
}

double calcPi_P1(int num_steps)
{
    double pi = 0.0;
    int inside_circle = 0;

    #pragma omp parallel reduction(+:inside_circle)
    {
        double x = 0;
        double y = 0;
        unsigned int seed = omp_get_thread_num() ^ time(NULL);
        #pragma omp for
        for(int i = 0; i < num_steps; i++)
        {
            x = ((double)rand_r(&seed)) / RAND_MAX;
            y = ((double)rand_r(&seed)) / RAND_MAX;
            double coords = (x*x) + (y*y);
            if(coords <= 1)
            {
                inside_circle++;
            }
        }
    }

    pi = 4.0 * ((double)inside_circle / (double)num_steps);
    return pi;
}


double calcPi_P2(int num_steps)
{
    double pi = 0.0;

    double area_sum = 0.0;
    double step = 2.0 / (double)num_steps;

    #pragma omp parallel reduction(+:area_sum)
    {
        #pragma omp for
        for(int i = 0; i < num_steps; i++)
        {
            double x = -1.0 + (i + 0.5) * step; 
            area_sum += sqrt(1.0 - (x*x));
            
        }
    }

    pi = 2.0 * step * area_sum;
    return pi;
}
