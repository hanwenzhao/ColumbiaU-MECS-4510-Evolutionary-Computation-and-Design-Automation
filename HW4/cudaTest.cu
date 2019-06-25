#include <stdio.h>

__global__ void vecAdd(void)
{
    printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
}

void vecAdd_wrapper() {

    // Execute the kernel
    vecAdd <<< 10, 10 >>> ();
    cudaDeviceSynchronize();
    printf("Hello World!\n");
}
