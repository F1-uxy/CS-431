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
void prefixwindow(int src[], int prefix[], uint32_t n);

typedef struct Node 
{
    int sum;
    struct Node *left;
    struct Node *right;
} Node;

void prefixtree(Node *node, int partialSum, int *prefixSum, int *idx);
Node* newNode(int sum);
Node* build(int *arr, int start, int end);
void prefixtree(Node *node, int partialSum, int *prefixSum, int *idx);

Node* newNode(int sum) 
{
    Node* node = (Node*)malloc(sizeof(Node));
    node->sum = sum;
    node->left = node->right = NULL;
    return node;
}


Node* build(int *arr, int start, int end) 
{
    if (start == end) {
        return newNode(arr[start]);
    }
    int mid = (start + end) / 2;
    Node* curr = newNode(0);

    #pragma omp task
    curr->left = build(arr, start, mid);

    #pragma omp task
    curr->right = build(arr, mid + 1, end);

    #pragma omp taskwait 
    curr->sum = curr->left->sum + curr->right->sum; 
    return curr;
}

void prefixtree(Node *node, int partialSum, int *prefixSum, int *idx) 
{
    if (!node->left && !node->right) { 
        #pragma omp critical
        {
            prefixSum[*idx] = partialSum + node->sum;
            (*idx)++;
        }
        return;
    }

    if(node->left != NULL)
    {
        #pragma omp task
        prefixtree(node->left, partialSum, prefixSum, idx);
    }
    
    if(node->right != NULL)
    {
        #pragma omp task
        prefixtree(node->right, partialSum + (node->left ? node->left->sum : 0), prefixSum, idx);
    }

    #pragma omp taskwait 
}

int main(int argc, char** argv)
{
    // get input values
    uint32_t arr_size = 100;
    uint32_t cpus = 0;
    if(argc > 1) {
        arr_size = atoi(argv[1]);
        cpus = atoi(argv[2]);
    } else {
        usage(argc, argv);
        printf("using %"PRIu32"\n", arr_size);
    }

    // Generate random array of size n
    int src[arr_size];
    int prefix[arr_size];
    memset(prefix, 0, sizeof(prefix));  
    srand(time(NULL));

    for(int r = 0; r < arr_size; r++)
    {
        int x = rand() % 10;
        src[r] = x;

    }

    // set up timer
    uint64_t start_t;
    uint64_t end_t;
    InitTSC();
    
    // create file for outputs
    FILE* results;
    char filename[32];
    sprintf(filename, "%d-%d.txt", cpus, arr_size);
    results = fopen(filename, "ab+");
    if (results == NULL) {
        printf("Cannot create file");
        return 1;
    }


    // calculate in serial 
    memcpy(prefix, src, arr_size * sizeof(int));
    start_t = ReadTSC();
    prefixSerial(src, prefix, arr_size);
    end_t = ReadTSC();

    printf("\nTime taken: %f \n", (ElapsedTime(end_t - start_t)));
    fprintf(results, "%g, ", ElapsedTime(end_t - start_t));
    
    // Calculate in parallel
    memcpy(prefix, src, arr_size * sizeof(int));
    start_t = ReadTSC();
    prefixwindow(src, prefix, arr_size);
    end_t = ReadTSC();

    printf("\nTime taken: %f \n", (ElapsedTime(end_t - start_t)));
    fprintf(results, "%g, ", ElapsedTime(end_t - start_t));



    Node* root = build(src, 0, arr_size - 1);
    int *sum = (int*)malloc(arr_size * sizeof(int));
    int idx = 0;
    start_t = ReadTSC();
    prefixtree(root, 0, sum, &idx);
    end_t = ReadTSC();

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

void prefixwindow(int src[], int prefix[], uint32_t n)
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
