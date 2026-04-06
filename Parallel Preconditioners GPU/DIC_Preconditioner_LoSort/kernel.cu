#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstring>
#include <vector>
using namespace std;

#define REPS 1000

__global__ void init
(
    const int nCells,
    const double* diagPtr,
    const double* upper,
    const int*  losortAddr,
    const int*  lowerAddr,
    const int* losortStartAddr,
    double* reciprocalDiagPtr
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nCells) {
        int end;
        int face;
        int index;
        int readCell;

        reciprocalDiagPtr[cell] = diagPtr[cell];

        __syncthreads();

        end = losortStartAddr[cell + 1];
        for (index = losortStartAddr[cell]; index < end; index++)
        {
            face = losortAddr[index];
            readCell = lowerAddr[index];
            reciprocalDiagPtr[cell] = upper[face] * upper[face] / reciprocalDiagPtr[readCell];
        }

        __syncthreads();

        reciprocalDiagPtr[cell] = 1.0 / reciprocalDiagPtr[cell];
    }
}

__global__ void precondition
(
    double* diagPtr,
    const int nCells,
    const double* upper,
    const int*  upperAddr,
    const int*  lowerAddr,
    const int*  losortAddr,
    const int* ownerStartAddr,
    const int* losortStartAddr,
    const double* initialResidual,
    const double* reciprocalDiagPtr
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nCells) {
        int end;
        int face;
        int index;
        int readCell;
        double tempValue;

        diagPtr[cell] = reciprocalDiagPtr[cell] * initialResidual[cell];

        __syncthreads();

        tempValue = 0.0;
        end = losortStartAddr[cell + 1];
        for (index = losortStartAddr[cell]; index < end; index++)
        {
            face = losortAddr[index];
            readCell = lowerAddr[index];
            tempValue += upper[face] * diagPtr[readCell];
        }
        diagPtr[cell] -= reciprocalDiagPtr[cell] * tempValue;

        __syncthreads();

        tempValue = 0.0;
        end = ownerStartAddr[cell + 1];
        for (face = ownerStartAddr[cell]; face < end; face++)
        {
            readCell = upperAddr[face];
            tempValue += upper[face] * diagPtr[readCell];
        }
        diagPtr[cell] -= reciprocalDiagPtr[cell] * tempValue;
    }
}

void destroy(
    int* sizes,
    double* diag,
    double* lower,
    double* upper,
    int* lowerAddr,
    int* upperAddr,
    int* ownerStartAddr,
    double* reciprocalDiagPtr,
    double* initialResidual,
    double* psi,
    int* losortAddr,
    int* losortStartAddr
) {
    delete[](sizes);
    delete[](diag);
    delete[](lower);
    delete[](upper);
    delete[](lowerAddr);
    delete[](upperAddr);
    delete[](ownerStartAddr);
    delete[](reciprocalDiagPtr);
    delete[](initialResidual);
    delete[](psi);
    delete[](losortAddr);
    delete[](losortStartAddr);
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
    int** ownerStartAddr,
    int** losortAddr,
    int** losortStartAddr
)
{
    FILE* binaryFile = fopen(fileName, "rb");

    *sizes = new int[10];
    fread(*sizes, sizeof(int), 10, binaryFile);
    *diag = new double[(*sizes)[0]];
    *lower = new double[(*sizes)[1]];
    *lowerAddr = new int[(*sizes)[2]];
    *upper = new double[(*sizes)[3]];
    *upperAddr = new int[(*sizes)[4]];
    *initialResidual = new double[(*sizes)[5]];
    *psi = new double[(*sizes)[6]];
    *ownerStartAddr = new int[(*sizes)[7]];
    *losortAddr = new int[(*sizes)[8]];
    *losortStartAddr = new int[(*sizes)[9]];

    fread(*diag, sizeof(double), (*sizes)[0], binaryFile);
    fread(*lower, sizeof(double), (*sizes)[1], binaryFile);
    fread(*lowerAddr, sizeof(int), (*sizes)[2], binaryFile);
    fread(*upper, sizeof(double), (*sizes)[3], binaryFile);
    fread(*upperAddr, sizeof(int), (*sizes)[4], binaryFile);
    fread(*initialResidual, sizeof(double), (*sizes)[5], binaryFile);
    fread(*psi, sizeof(double), (*sizes)[6], binaryFile);
    fread(*ownerStartAddr, sizeof(int), (*sizes)[7], binaryFile);
    fread(*losortAddr, sizeof(int), (*sizes)[8], binaryFile);
    fread(*losortStartAddr, sizeof(int), (*sizes)[9], binaryFile);

    fclose(binaryFile);
}

int main(int argc, char* argv[]) {
    // CPU variables
    int* sizes;
    double* diag;
    double* diagCopy;
    double* lower;
    double* upper;
    int* lowerAddr;
    int* upperAddr;
    int* losortAddr;
    int* ownerStartAddr;
    int* losortStartAddr;
    double* reciprocalDiagPtr;
    double* reciprocalDiagPtrCopy;
    double* initialResidual;
    double* psi;
    // GPU variables
    int* dev_upperAddr;
    int* dev_lowerAddr;
    int* dev_losortAddr;
    int* dev_ownerStartAddr;
    int* dev_losortStartAddr;
    double* dev_diag;
    double* dev_upper;
    double* dev_initialResidual;
    double* dev_reciprocalDiagPtr;
    // CUDA time measuring variables
    cudaEvent_t start;
    cudaEvent_t stop;
    float elapsedTime;

    cudaSetDevice(0);

    //Read matrix values
    binaryRead("./matrixBeforeLong.bin", &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi, &ownerStartAddr, &losortAddr, &losortStartAddr);

    // Reorder lowerAddr to avoid double memory acces with losortAddr
    int* newLowerAddr = new int[sizes[1]];
    for (int cell = 0; cell < sizes[0]; cell++)
    {
        int end = losortStartAddr[cell + 1];
        for (int index = losortStartAddr[cell]; index < end; index++)
        {
            int face = losortAddr[index];
            newLowerAddr[index] = lowerAddr[face];
        }
    }
    delete[](lowerAddr);
    lowerAddr = newLowerAddr;

    // Create time events
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Reserve GPU mem
    cudaMalloc((void**)&dev_upperAddr, sizes[1] * sizeof(int));
    cudaMalloc((void**)&dev_lowerAddr, sizes[1] * sizeof(int));
    cudaMalloc((void**)&dev_diag, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_losortAddr, sizes[1] * sizeof(int));
    cudaMalloc((void**)&dev_ownerStartAddr, (sizes[0] + 1) * sizeof(int));
    cudaMalloc((void**)&dev_losortStartAddr, (sizes[0] + 1) * sizeof(int));
    cudaMalloc((void**)&dev_upper, sizes[1] * sizeof(double));
    cudaMalloc((void**)&dev_initialResidual, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_reciprocalDiagPtr, sizes[0] * sizeof(double));

    // Send data to GPU
    cudaMemcpy(dev_upperAddr, upperAddr, sizes[1] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lowerAddr, lowerAddr, sizes[1] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_diag, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_losortAddr, losortAddr, sizes[1] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ownerStartAddr, ownerStartAddr, (sizes[0] + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_losortStartAddr, losortStartAddr, (sizes[0] + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_upper, upper, sizes[1] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_initialResidual, initialResidual, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);

    // Checking vectors config
    diagCopy = new double[sizes[0]];

    reciprocalDiagPtr = new double[sizes[0]];
    reciprocalDiagPtrCopy = new double[sizes[0]];

    for (int k = 8; k <= 1024; k *= 2)
    {
        int threadsPerBlock = k;
        int blocksPerGrid = (sizes[0] + threadsPerBlock - 1) / threadsPerBlock;

        cout << "DIC Preconditioner CUDA (" << k << " Num Threads):" << endl;

        for (int j = 0; j < 2; j++)
        {
            cudaMemcpy(dev_diag, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();

            // Time init
            cudaEventRecord(start, 0);
            for (int i = 0; i < REPS; i++) {
                init << <blocksPerGrid, threadsPerBlock >> >
                        (sizes[0], dev_diag, dev_upper, dev_losortAddr, dev_lowerAddr, dev_losortStartAddr, dev_reciprocalDiagPtr);                
                cudaDeviceSynchronize();
            }
            cudaEventRecord(stop, 0);
            cudaEventSynchronize(stop);

            // Get time in seconds
            cudaEventElapsedTime(&elapsedTime, start, stop);
            cout << (elapsedTime / REPS) / 1000 << endl;

            // Time precondition
            cudaEventRecord(start, 0);
            for (int i = 0; i < REPS; i++) {
                precondition << <blocksPerGrid, threadsPerBlock >> > (dev_diag, sizes[0], dev_upper, dev_upperAddr, dev_lowerAddr, 
                                 dev_losortAddr, dev_ownerStartAddr, dev_losortStartAddr, dev_initialResidual, dev_reciprocalDiagPtr);
                cudaDeviceSynchronize();
            }
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
    cout << "DIC Preconditioner CUDA (Sequential):" << endl;
    clock_t t0, t1;
    for (int j = 0; j < 2; j++)
    {
        // Fill diagCopy with original values
        for (int cell = 0; cell < sizes[0]; cell++)
        {
            diagCopy[cell] = diag[cell];
        }

        t0 = clock();
        for (int i = 0; i < 100; i++)
        {
            for (int cell = 0; cell < sizes[0]; cell++)
            {
                reciprocalDiagPtrCopy[cell] = diagCopy[cell];
            }

            for (int face = 0; face < sizes[1]; face++)
            {
                reciprocalDiagPtrCopy[upperAddr[face]] -= upper[face] * upper[face] / reciprocalDiagPtrCopy[lowerAddr[face]];
            }

            for (int cell = 0; cell < sizes[0]; cell++)
            {
                reciprocalDiagPtrCopy[cell] = 1.0 / reciprocalDiagPtrCopy[cell];
            }
        }
        t1 = clock();
        cout << (double(t1 - t0) / CLOCKS_PER_SEC) / 100 << endl;

        // Precondition
        t0 = clock();
        int nFacesM1 = sizes[1] - 1;
        for (int i = 0; i < 100; i++)
        {
            for (int cell = 0; cell < sizes[0]; cell++)
            {
                diagCopy[cell] = reciprocalDiagPtrCopy[cell] * initialResidual[cell];
            }

            for (int face = 0; face < sizes[1]; face++)
            {
                diagCopy[upperAddr[face]] -= reciprocalDiagPtrCopy[upperAddr[face]] * upper[face] * diagCopy[lowerAddr[face]];
            }

            for (int face = nFacesM1; face >= 0; face--)
            {
                diagCopy[lowerAddr[face]] -= reciprocalDiagPtrCopy[lowerAddr[face]] * upper[face] * diagCopy[upperAddr[face]];
            }
        }
        t1 = clock();
        cout << (double(t1 - t0) / CLOCKS_PER_SEC) / 100 << endl;
    }

    // Return to CPU result data to compare
    cudaMemcpy(diag, dev_diag, sizes[0] * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(reciprocalDiagPtr, dev_reciprocalDiagPtr, sizes[0] * sizeof(double), cudaMemcpyDeviceToHost);

    // Free GPU mem
    cudaFree(dev_upperAddr);
    cudaFree(dev_lowerAddr);
    cudaFree(dev_diag);
    cudaFree(dev_losortAddr);
    cudaFree(dev_losortStartAddr);
    cudaFree(dev_upper);
    cudaFree(dev_initialResidual);
    cudaFree(dev_reciprocalDiagPtr);
    cudaFree(dev_ownerStartAddr);


    // Checking
    /*for (int j = 0; j < sizes[0]; j++)
    {

        if (fabs(reciprocalDiagPtr[j] - reciprocalDiagPtrCopy[j]) > 1e-20)
        {
            cout << "Index: " << j << " error in reciprocalDiagPtr == " << reciprocalDiagPtr[j] << " and reciprocalDiagPtrCopy == " << reciprocalDiagPtrCopy[j] << endl;
            break;
        }
        if (fabs(diag[j] - diagCopy[j]) > 1e-20)
        {
            cout << "Index: " << j << " error in diag == " << diag[j] << " and diagCopy == " << diagCopy[j] << endl;
            break;
        }
    }*/

    delete[](diagCopy);
    delete[](reciprocalDiagPtrCopy);

    //Destroy matrix
    destroy(sizes, diag, lower, upper, lowerAddr, upperAddr, ownerStartAddr, reciprocalDiagPtr, initialResidual, psi, losortAddr, losortStartAddr);

    return 0;
}