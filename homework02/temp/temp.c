#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

void blelloch_scan(int *arr, int *output, int n) {
    int *temp = (int*) malloc(n * sizeof(int));
    int *arr2 = (int*) malloc(n * sizeof(int));

    
    // Copy the original array to a temporary array for the final adjustment
    for (int i = 0; i < n; i++) {
        temp[i] = arr[i];
        //printf("%d vs %d, ", arr[i], temp[i]);
    }

    // Step 1: Upsweep phase (reduce)
    int d, i;
    for (d = 1; d < n; d *= 2) {
        for (i = 0; i < n; i += 2 * d) {
            if (i + d < n) {
                arr[i + 2 * d - 1] += arr[i + d - 1];
            }
        }
    }
    
    // Set the last element to 0 to start the down-sweep
    arr[n - 1] = 0;

    // Step 2: Downsweep phase
    for (d = n / 2; d >= 1; d /= 2) {
        for (i = 0; i < n; i += 2 * d) {
            if (i + d < n) {
                int tmp = arr[i + d - 1];
                arr[i + d - 1] = arr[i + 2 * d - 1];
                arr[i + 2 * d - 1] += tmp;
            }
        }
    }

    // Copy the result to the output array and add the original values
    for (i = 0; i < n; i++) {
        //printf("%d, ", arr2[i]);
        output[i] = arr[i] + temp[i];
    }

    free(temp);
}


int main() {

    int src[8];
    srand(time(NULL));

    for(int r = 0; r < 8; r++)
    {
        int x = rand() % 10;
        src[r] = x;
        printf("%d, ", x);

    }

    int arr[] = {3, 1, 7, 0, 4, 1, 6, 3};
    int arr2[] = {3, 1, 7, 0, 4, 1, 6, 3};

    int n = sizeof(src) / sizeof(src[0]);
    int *output = (int*) malloc(n * sizeof(int));

    blelloch_scan(src, output, n);

    printf("Prefix Sum Result:\n");

    
    for (int i = 0; i < n; i++) {
        printf("%d ", output[i]);
    }
    printf("\n");

    free(output);
    return 0;
}
