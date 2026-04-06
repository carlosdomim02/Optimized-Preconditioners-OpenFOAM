#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
using namespace std;

#define REPS 1000000

__global__ void init
(
    const int nCells,
    const double* diagPtr,
    double* reciprocalDiagPtr
)
{
    for (int i = 0; i < REPS; i++) {
        int cell = blockDim.x * blockIdx.x + threadIdx.x;

        if (cell < nCells) {
            reciprocalDiagPtr[cell] = 1.0 / diagPtr[cell];
        }
    }
}

__global__ void precondition
(
    double* diag,
    const int nCells,
    const double* initialResidual,
    const double* reciprocalDiagPtr
)
{
    for (int i = 0; i < REPS; i++) {
        int cell = blockDim.x * blockIdx.x + threadIdx.x;
            
        if (cell < nCells) {
            diag[cell] = reciprocalDiagPtr[cell] * initialResidual[cell];
        }
    }
}

void initSeq
(
    int nCells,
    double* diagPtr,
    double* reciprocalDiagPtr
)
{
    for (int cell = 0; cell < nCells; cell++)
    {
        reciprocalDiagPtr[cell] = 1.0 / diagPtr[cell];
    }
}

void preconditionSeq
(
    double* diag,
    const int nCells,
    const double* initialResidual,
    const double* reciprocalDiagPtr
)
{
    for (int cell = 0; cell < nCells; cell++)
    {
        diag[cell] = reciprocalDiagPtr[cell] * initialResidual[cell];
    }
}

void destroy
(
    const int* sizes,
    const double* psi,
    const double* diag,
    const double* upper,
    const double* lower,
    const int* upperAddr,
    const int* lowerAddr,
    const double* initialResidual,
    const int* ownerStartAddr,
    const double* reciprocalDiagPtr
)
{
    delete[](psi);
    delete[](diag);
    delete[](sizes);
    delete[](upper);
    delete[](lower);
    delete[](upperAddr);
    delete[](lowerAddr);
    delete[](ownerStartAddr);
    delete[](initialResidual);
}

void binaryRead
(
    char* fileName,
    int** sizes,
    double** diag,
    int** upperAddr,
    int** lowerAddr,
    double** upper,
    double** lower,
    double** initialResidual,
    double** psi,
    int** ownerStartAddr
)
{
    FILE* binaryFile = fopen(fileName, "rb");

    *sizes = new int[8];
    fread(*sizes, sizeof(int), 8, binaryFile);
    *diag = new double[(*sizes)[0]];
    *lower = new double[(*sizes)[1]];
    *lowerAddr = new int[(*sizes)[2]];
    *upper = new double[(*sizes)[3]];
    *upperAddr = new int[(*sizes)[4]];
    *initialResidual = new double[(*sizes)[5]];
    *psi = new double[(*sizes)[6]];
    *ownerStartAddr = new int[(*sizes)[7]];

    fread(*diag, sizeof(double), (*sizes)[0], binaryFile);
    fread(*lower, sizeof(double), (*sizes)[1], binaryFile);
    fread(*lowerAddr, sizeof(int), (*sizes)[2], binaryFile);
    fread(*upper, sizeof(double), (*sizes)[3], binaryFile);
    fread(*upperAddr, sizeof(int), (*sizes)[4], binaryFile);
    fread(*initialResidual, sizeof(double), (*sizes)[5], binaryFile);
    fread(*psi, sizeof(double), (*sizes)[6], binaryFile);
    fread(*ownerStartAddr, sizeof(int), (*sizes)[7], binaryFile);

    fclose(binaryFile);
}

int main() {
    // CPU variables
    int* sizes;
    double* psi;
    double* diag;
    double* diagCopy;
    double* lower;
    double* upper;
    int* lowerAddr;
    int* upperAddr;
    int* ownerStartAddr;
    double* initialResidual;
    double* reciprocalDiagPtr;
    double* reciprocalDiagPtrCopy;
    // GPU variables
    double* dev_diag;
    double* dev_initialResidual;
    double* dev_reciprocalDiagPtr;
    cudaDeviceProp prop;
    // Time measure variables
    cudaEvent_t start;
    cudaEvent_t stop;
    float elapsedTime;

    cudaSetDevice(1);

    //Read matrix values
    binaryRead("./matrixBefore.bin", &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi, &ownerStartAddr);

    // Create time events
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Reserve GPU mem
    cudaMalloc((void**)&dev_diag, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_initialResidual, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_reciprocalDiagPtr, sizes[0] * sizeof(double));

    // Send data to GPU
    cudaMemcpy(dev_diag, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_initialResidual, initialResidual, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);

    // Checking vectors config
    diagCopy = new double[sizes[0]];

    reciprocalDiagPtr = new double[sizes[0]];
    reciprocalDiagPtrCopy = new double[sizes[0]];

    for (int k = 8; k <= 1024; k *= 2)
    {
        int threadsPerBlock = k;
        int blocksPerGrid = (sizes[0] + threadsPerBlock - 1) / threadsPerBlock;

        cout << "Diag Preconditioner CUDA (" << k << " Num Threads):" << endl;

        for (int j = 0; j < 2; j++)
        {
            cudaMemcpy(dev_diag, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();

            // Time init
            cudaEventRecord(start, 0);
            init << <blocksPerGrid, threadsPerBlock >> > (sizes[0], dev_diag, dev_reciprocalDiagPtr);
            cudaDeviceSynchronize();
            cudaDeviceSynchronize();
            cudaEventRecord(stop, 0);
            cudaEventSynchronize(stop);

            // Get time in seconds
            cudaEventElapsedTime(&elapsedTime, start, stop);
            cout << (elapsedTime / REPS) / 1000 << endl;

            // Time precondition
            cudaEventRecord(start, 0);
            precondition << <blocksPerGrid, threadsPerBlock >> > (dev_diag, sizes[0], dev_initialResidual, dev_reciprocalDiagPtr);
            cudaDeviceSynchronize();
            cudaEventRecord(stop, 0);
            cudaEventSynchronize(stop);

            // Get time in seconds
            cudaEventElapsedTime(&elapsedTime, start, stop);
            cout << (elapsedTime / REPS) / 1000 << endl;
        }
    }

    // Free streams and events
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // Check results   
    // Init 
    cout << "Diag Preconditioner CUDA (Secuencial):" << endl;
    clock_t t0, t1;
    for (int j = 0; j < 2; j++)
    {
        // Fill diagCopy with original values
        for (int cell = 0; cell < sizes[0]; cell++)
        {
            diagCopy[cell] = diag[cell];
        }

        t0 = clock();
        for (int i = 0; i < 10000; i++)
        {
            initSeq(sizes[0], diagCopy, reciprocalDiagPtrCopy);
        }
        t1 = clock();
        cout << (double(t1 - t0) / CLOCKS_PER_SEC) / 10000 << endl;

        // Precondition
        t0 = clock();
        for (int i = 0; i < 10000; i++)
        {
            preconditionSeq(diagCopy, sizes[0], initialResidual, reciprocalDiagPtrCopy);
        }
        t1 = clock();
        cout << (double(t1 - t0) / CLOCKS_PER_SEC) / 10000 << endl;
    }

    // Return to CPU result data to compare
    cudaMemcpy(diag, dev_diag, sizes[0] * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(reciprocalDiagPtr, dev_reciprocalDiagPtr, sizes[0] * sizeof(double), cudaMemcpyDeviceToHost);

    // Free GPU mem
    cudaFree(dev_diag);
    cudaFree(dev_initialResidual);
    cudaFree(dev_reciprocalDiagPtr);

    // Checking
    for(int j = 0; j < sizes[0]; j++)
    {

        if(fabs(reciprocalDiagPtr[j] - reciprocalDiagPtrCopy[j]) > 1e-20)
        {
            cout << "Index: " << j << " error in reciprocalDiagPtr == " << reciprocalDiagPtr[j] << " and reciprocalDiagPtrCopy == " << reciprocalDiagPtrCopy[j] << endl;
            break;
        }
        if(fabs(diag[j] - diagCopy[j]) > 1e-20)
        {
            cout << "Index: " << j << " error in diag == " << diag[j] << " and diagCopy == " << diagCopy[j] << endl;
            break;
        }
    }

    delete[](diagCopy);
    delete[](reciprocalDiagPtrCopy);

    // Destroy matrix
    destroy(sizes, psi, diag, upper, lower, upperAddr, lowerAddr, initialResidual, ownerStartAddr, reciprocalDiagPtr);

    return 0;
}