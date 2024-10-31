#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

typedef struct TreeNode {
    int sum;
    struct TreeNode *left, *right;
} TreeNode;

TreeNode* createNode(int sum) {
    TreeNode* node = (TreeNode*)malloc(sizeof(TreeNode));
    node->sum = sum;
    node->left = node->right = NULL;
    return node;
}

TreeNode* buildTree(int *arr, int start, int end) {
    if (start == end) {
        return createNode(arr[start]);
    }
    int mid = (start + end) / 2;
    TreeNode* root = createNode(0);

    #pragma omp task shared(root)
    root->left = buildTree(arr, start, mid);

    #pragma omp task shared(root)
    root->right = buildTree(arr, mid + 1, end);

    #pragma omp taskwait 
    root->sum = root->left->sum + root->right->sum; 
    return root;
}

void calculatePrefixSums(TreeNode *node, int partialSum, int *prefixSum, int *index) {
    if (!node->left && !node->right) { 
        #pragma omp critical
        {
            prefixSum[*index] = partialSum + node->sum;
            (*index)++;
        }
        return;
    }

    #pragma omp task if(node->left)  
    calculatePrefixSums(node->left, partialSum, prefixSum, index);

    #pragma omp task if(node->right) 
    calculatePrefixSums(node->right, partialSum + (node->left ? node->left->sum : 0), prefixSum, index);

    #pragma omp taskwait 
}

void freeTree(TreeNode *node) {
    if (node == NULL) return;
    freeTree(node->left);
    freeTree(node->right);
    free(node);
}

int main() {
    int arr[] = {1, 2, 3, 4, 5, 6, 7, 8};
    int n = sizeof(arr) / sizeof(arr[0]);

    // Step 1: Build the binary tree
    TreeNode* root;
    #pragma omp parallel
    {
        #pragma omp single
        root = buildTree(arr, 0, n - 1);
    }

    // Step 2: Calculate the prefix sums
    int *prefixSum = (int*)malloc(n * sizeof(int));
    int index = 0;
    #pragma omp parallel
    {
        #pragma omp single
        calculatePrefixSums(root, 0, prefixSum, &index);
    }

    // Print the prefix sums
    printf("Prefix Sums: ");
    for (int i = 0; i < n; i++) {
        printf("%d ", prefixSum[i]);
    }
    printf("\n");

    // Free the tree and prefix sum array
    freeTree(root);
    free(prefixSum);

    return 0;
}
