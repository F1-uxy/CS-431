#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "common.h"
#include <inttypes.h>
#include <time.h>
#include <string.h>


void usage(int argc, char** argv);
void prefixSerial(int src[], int prefix[], uint32_t n);
void prefixlogn(int src[], int prefix[], uint32_t n);
void prefixn(int src[], int prefix[], uint32_t n);


int main(int argc, char** argv)
{
    // get input values
    uint32_t arr_size = 100;
    if(argc > 1) {
        arr_size = atoi(argv[1]);
    } else {
        usage(argc, argv);
        printf("using %"PRIu32"\n", arr_size);
    }

    // Generate random array of size n
    int src[arr_size];
    int src2[arr_size];
    int prefix[arr_size];
    int prefix2[arr_size];
    memset(prefix, 0, sizeof(prefix));  
    srand(time(NULL));

    for(int r = 0; r < arr_size; r++)
    {
        int x = rand() % 10;
        src[r] = x;
        src2[r] = x;

    }
    printf("\n");

    // set up timer
    uint64_t start_t;
    uint64_t end_t;
    InitTSC();
    
    // create file for outputs
    FILE* results;
    char filename[32];
    sprintf(filename, "%d.txt", arr_size);
    results = fopen(filename, "ab+");
    if (results == NULL) {
        printf("Cannot create file");
        return 1;
    }

    for(int i = 0; i < arr_size; i++)
    {
        printf("%d, ", src[i]);
    }
    printf("\n");


    // calculate in serial 
    memcpy(prefix, src, arr_size * sizeof(int));
    start_t = ReadTSC();
    prefixSerial(src, prefix, arr_size);
    end_t = ReadTSC();
    for(int i = 0; i < arr_size; i++)
    {
        printf("%d, ", prefix[i]);
    }
    
    printf("\nTime taken: %f \n", (ElapsedTime(end_t - start_t)));
    fprintf(results, "%g, ", ElapsedTime(end_t - start_t));
    


    memcpy(prefix, src, arr_size * sizeof(int));
    start_t = ReadTSC();
    prefixlogn(src, prefix, arr_size);
    end_t = ReadTSC();
    for(int i = 0; i < arr_size; i++)
    {
        printf("%d, ", prefix[i]);
    }
    printf("\nTime taken: %f \n", (ElapsedTime(end_t - start_t)));
    fprintf(results, "%g, ", ElapsedTime(end_t - start_t));



    int *output = (int*)malloc(arr_size * sizeof(int));
    start_t = ReadTSC();
    prefixn(src, output, arr_size);
    end_t = ReadTSC();

    for (int i = 0; i < arr_size; i++) {
        printf("%d, ", output[i]);
    }
    printf("\nTime taken: %f \n", ElapsedTime(end_t - start_t));
    fprintf(results, "%g \n", ElapsedTime(end_t - start_t));

    
    
    return 0;
}


void usage(int argc, char** argv)
{
    fprintf(stdout, "usage: %s <# arr size>\n", argv[0]);
}

void prefixSerial(int src[], int prefix[], uint32_t n)
{
    for (uint32_t i = 0; i < n; i++)
    {
        prefix[i] = 0;
        for (uint32_t j = 0; j <= i; j++)
        {
            prefix[i] += src[j];
        }
    }
}

void prefixlogn(int src[], int prefix[], uint32_t n)
{
    int temp[n];

    for (int i = 1; i < n; i *= 2)
    {
        memcpy(temp, prefix, n * sizeof(int)); // Remove race condition to prefix by copying array

        #pragma omp parallel for
        for (int j = i; j < n; j++)
        {
            prefix[j] = temp[j] + temp[j - i];
        }
    }
}

// Using Belloch Scan algorithm to simulate O(n) prefix sum
void prefixn(int *arr, int *output, uint32_t n)
{
    int *temp = (int*) malloc(n * sizeof(int));
    int *arr2 = (int*) malloc(n * sizeof(int));

    
    for (int i = 0; i < n; i++) {
        temp[i] = arr[i];
    }

    int d, i;
    for (d = 1; d < n; d *= 2) {
        #pragma omp parallel for
        for (i = 0; i < n; i += 2 * d) {
            if (i + d < n) {
                arr[i + 2 * d - 1] += arr[i + d - 1];
            }
        }
    }
    
    arr[n - 1] = 0;

    for (d = n / 2; d >= 1; d /= 2) {
        #pragma omp parallel for
        for (i = 0; i < n; i += 2 * d) {
            if (i + d < n) {
                int tmp = arr[i + d - 1];
                arr[i + d - 1] = arr[i + 2 * d - 1];
                arr[i + 2 * d - 1] += tmp;
            }
        }
    }

    for (i = 0; i < n; i++) {
        output[i] = arr[i] + temp[i];
    }

    free(temp);
}
