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

// Waiting Version
__global__ void init1
(
    const int nCells,
    const double* diagPtr,
    double* reciprocalDiagPtr,
    int* dependenciesCopy
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nCells) {
        reciprocalDiagPtr[cell] = diagPtr[cell];
        dependenciesCopy[cell] = 0;
    }
}

__device__ double atomicAddDouble(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);

    return __longlong_as_double(old);
}

__global__ void init2(double* reciprocalDiagPtr, const double* upper, const int* lowerAddr, const int* upperAddr, int* dependenciesCopy, const int* dependencies, int nFaces) {
    int face = blockIdx.x * blockDim.x + threadIdx.x;
    double expr;

    if (face < nFaces) {
        // Espera activa mientras se cumplen las dependencias
        while (dependencies[lowerAddr[face]] > dependenciesCopy[lowerAddr[face]]) {
            // Puede incluir algún tipo de espera activa o pausa si es necesario.
        }

        expr = upper[face] * upper[face] / reciprocalDiagPtr[lowerAddr[face]];

        // Actualiza atomícamente el valor
        atomicAddDouble(&(reciprocalDiagPtr[upperAddr[face]]), -expr);
        //atomicAdd(&(reciprocalDiagPtr[upperAddr[face]]), -expr);
        atomicAdd(&(dependenciesCopy[upperAddr[face]]), 1);
    }
}

__global__ void init3
(
    const int nCells,
    double* reciprocalDiagPtr
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nCells) {
        reciprocalDiagPtr[cell] = 1.0 / reciprocalDiagPtr[cell];
    }
}

__global__ void precondition1
(
    double* diag,
    const int nCells,
    const double* initialResidual,
    const double* reciprocalDiagPtr,
    int* dependenciesCopy,
    int* dependenciesReverseCopy
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nCells) {
        diag[cell] = reciprocalDiagPtr[cell] * initialResidual[cell];
        dependenciesCopy = 0;
        dependenciesReverseCopy = 0;
    }
}

__global__ void precondition2(double* reciprocalDiagPtr, double* diagPtr, const double* upper, const int* lowerAddr, const int* upperAddr, int* dependenciesCopy, const int* dependencies, int nFaces) {
    int face = blockIdx.x * blockDim.x + threadIdx.x;
    double expr;

    if (face < nFaces) {
        // Espera activa mientras se cumplen las dependencias
        while (dependencies[lowerAddr[face]] > dependenciesCopy[lowerAddr[face]]) {
            // Puede incluir algún tipo de espera activa o pausa si es necesario.
        }

        expr = reciprocalDiagPtr[upperAddr[face]] * upper[face] * diagPtr[lowerAddr[face]];
        // Actualiza atomícamente el valor
        atomicAddDouble(&(diagPtr[upperAddr[face]]), -expr);
        //atomicAdd(&(reciprocalDiagPtr[upperAddr[face]]), -expr);
        atomicAdd(&(dependenciesCopy[upperAddr[face]]), 1);
    }
}

__global__ void precondition3(double* reciprocalDiagPtr, double* diagPtr, const double* upper, const int* lowerAddr, const int* upperAddr, int* dependenciesCopy, const int* dependencies, int nFaces) {
    int face = blockIdx.x * blockDim.x + threadIdx.x;
    double expr;

    if (face < nFaces) {
        // Espera activa mientras se cumplen las dependencias
        while (dependencies[upperAddr[face]] > dependenciesCopy[upperAddr[face]]) {
            // Puede incluir algún tipo de espera activa o pausa si es necesario.
        }

        expr = reciprocalDiagPtr[lowerAddr[face]] * upper[face] * diagPtr[upperAddr[face]];

        // Actualiza atomícamente el valor
        atomicAddDouble(&(diagPtr[lowerAddr[face]]), -expr);
        //atomicAdd(&(reciprocalDiagPtr[upperAddr[face]]), -expr);
        atomicAdd(&(dependenciesCopy[lowerAddr[face]]), 1);
    }
}

void init
(
    int nCells,
    int nFaces,
    double* upper,
    int* upperAddr,
    int* lowerAddr,
    double* diagPtr,
    double** reciprocalDiagPtr
)
{
    (*reciprocalDiagPtr) = new double[nCells];

    //Initialization of reciprocalDiag
    for (int cell = 0; cell < nCells; cell++)
    {
        (*reciprocalDiagPtr)[cell] = diagPtr[cell];
    }

    // Calculate the DIC diagonal
    for (int face = 0; face < nFaces; face++)
    {
        (*reciprocalDiagPtr)[upperAddr[face]] -= upper[face] * upper[face] / (*reciprocalDiagPtr)[lowerAddr[face]];
    }

    // Calculate the reciprocal of the preconditioned diagonal
    for (int cell = 0; cell < nCells; cell++)
    {
        (*reciprocalDiagPtr)[cell] = 1.0 / (*reciprocalDiagPtr)[cell];
    }
}

void precondition
(
    double* diagPtr,
    const int nCells,
    const int nFaces,
    const double* upper,
    const int* upperAddr,
    const int* lowerAddr,
    const double* initialResidual,
    const double* reciprocalDiagPtr
)
{
    int nFacesM1 = nFaces - 1;

    for (int cell = 0; cell < nCells; cell++)
    {
        diagPtr[cell] = reciprocalDiagPtr[cell] * initialResidual[cell];
    }

    for (int face = 0; face < nFaces; face++)
    {
        diagPtr[upperAddr[face]] -= reciprocalDiagPtr[upperAddr[face]] * upper[face] * diagPtr[lowerAddr[face]];
    }

    for (int face = nFacesM1; face >= 0; face--)
    {
        diagPtr[lowerAddr[face]] -= reciprocalDiagPtr[lowerAddr[face]] * upper[face] * diagPtr[upperAddr[face]];
    }
}

void destroy(
    const int* sizes,
    double* diag,
    double* lower,
    double* upper,
    int* lowerAddr,
    int* upperAddr,
    int* ownerStartAddr,
    double* initialResidual,
    double* psi,
    int* dependencies,
    int* dependenciesReverse)
{
    delete[] sizes;
    delete[] diag;
    delete[] lower;
    delete[] upper;
    delete[] lowerAddr;
    delete[] upperAddr;
    delete[] ownerStartAddr;
    delete[] initialResidual;
    delete[] psi;
    delete[] dependencies;
    delete[] dependenciesReverse;
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

void binaryPrint
(
    const char* fileName,
    const int* sizes,
    const double* diag,
    const int* upperAddr,
    const int* lowerAddr,
    const double* upper,
    const double* lower,
    const double* initialResidual,
    const double* psi,
    const int* ownerStartAddr
)
{
    FILE* binaryFile = fopen(fileName, "wb");

    fwrite(sizes, sizeof(int), 8, binaryFile);
    fwrite(diag, sizeof(double), sizes[0], binaryFile);
    fwrite(lower, sizeof(double), sizes[1], binaryFile);
    fwrite(lowerAddr, sizeof(int), sizes[2], binaryFile);
    fwrite(upper, sizeof(double), sizes[3], binaryFile);
    fwrite(upperAddr, sizeof(int), sizes[4], binaryFile);
    fwrite(initialResidual, sizeof(double), sizes[5], binaryFile);
    fwrite(psi, sizeof(double), sizes[6], binaryFile);
    fwrite(ownerStartAddr, sizeof(int), sizes[7], binaryFile);

    fclose(binaryFile);
}

void orderDependencies
(
    const int nCells,
    const int nFaces,
    const double* upper,
    const int* upperAddr,
    const int* lowerAddr,
    const int* ownerStartAddr,
    int** dependencies
)
{
    (*dependencies) = new int[nCells];
    memset((*dependencies), 0, nCells * sizeof(int));

    for (int i = 0; i < nFaces; i++) {
        if (ownerStartAddr[upperAddr[i] + 1] - ownerStartAddr[upperAddr[i]] > 0) (*dependencies)[upperAddr[i]]++;
        else (*dependencies)[upperAddr[i]] = -1;
    }
}

void orderDependenciesReverse
(
    const int nCells,
    const int nFaces,
    const double* upper,
    const int* upperAddr,
    const int* lowerAddr,
    const int* ownerStartAddr,
    int** dependencies
)
{
    (*dependencies) = new int[nCells];
    memset((*dependencies), 0, nCells * sizeof(int));

    for (int i = 0; i < nFaces; i++) {
        if (ownerStartAddr[lowerAddr[i] + 1] - ownerStartAddr[lowerAddr[i]] > 0) (*dependencies)[lowerAddr[i]]++;
        else (*dependencies)[lowerAddr[i]] = -1;
    }
}

int main() {
    // CPU variables
    int* sizes;
    double* diag;
    double* diagCopy;
    double* lower;
    double* upper;
    double* reciprocalDiagPtr;
    int* lowerAddr;
    int* upperAddr;
    int* ownerStartAddr;
    int* dependencies;
    int* dependenciesReverse;
    int* dependenciesReverseCopy;
    double* initialResidual;
    double* psi;
    // GPU variables
    int* dev_upperAddr;
    int* dev_lowerAddr;
    int* dev_dependencies;
    int* dev_dependenciesCopy;
    int* dev_dependenciesReverse;
    int* dev_dependenciesReverseCopy;
    double* dev_diag;
    double* dev_diagCopy;
    double* dev_upper;
    double* dev_initialResidual;
    double* dev_reciprocalDiagPtr;
    // Special cuda variables
    cudaEvent_t start;
    cudaEvent_t stop;
    float elapsedTime;
    // File variables
    ofstream outputFile;
    outputFile.open("output.txt");

    cudaSetDevice(1);

    //Read matrix values
    binaryRead("./matrixBefore.bin", &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi, &ownerStartAddr);

    //Ordder Dependencies
    orderDependencies(sizes[0], sizes[1], upper, upperAddr, lowerAddr, ownerStartAddr, &dependencies);
    orderDependenciesReverse(sizes[0], sizes[1], upper, upperAddr, lowerAddr, ownerStartAddr, &dependenciesReverse);

    // Create events and streams
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Reserve GPU mem
    cudaMalloc((void**)&dev_upperAddr, sizes[3] * sizeof(int));
    cudaMalloc((void**)&dev_lowerAddr, sizes[3] * sizeof(int));
    cudaMalloc((void**)&dev_upper, sizes[3] * sizeof(double));
    cudaMalloc((void**)&dev_initialResidual, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_diag, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_diagCopy, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_reciprocalDiagPtr, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_dependencies, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_dependenciesCopy, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_dependenciesReverse, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_dependenciesReverseCopy, sizes[0] * sizeof(double));

    // Send data to GPU
    cudaMemcpy(dev_upperAddr, upperAddr, sizes[3] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lowerAddr, lowerAddr, sizes[3] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_upper, upper, sizes[3] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_initialResidual, initialResidual, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_diag, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_diagCopy, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_dependencies, dependencies, sizes[0] * sizeof(int), cudaMemcpyHostToDevice);

    // Secuencial option
    /*cudaDeviceSynchronize();
    clock_t t0, t1;
    diagCopy = new double[sizes[0]];
    t0 = clock();
    for (int i = 0; i < 10000; i++) {
        init(sizes[0], sizes[3], upper, upperAddr, lowerAddr, diag, &reciprocalDiagPtr);
    }
    t1 = clock();
    printf("Tiempo init secuencial: %f s\n", (double(t1 - t0) / CLOCKS_PER_SEC) / 10000);

    t0 = clock();
    for (int i = 0; i < 10000; i++) {
        precondition(diag, sizes[0], sizes[3], upper, upperAddr, lowerAddr, initialResidual, reciprocalDiagPtr);
    }
    t1 = clock();
    printf("Tiempo precondition secuencial: %f s\n", (double(t1 - t0) / CLOCKS_PER_SEC) / 10000);*/

    // Define number of blocks and threads
    int threadsPerBlock = 64;
    int blocksPerGrid = (sizes[0] + threadsPerBlock - 1) / threadsPerBlock;
    int blocksPerKernel;

    //Init
    cudaDeviceSynchronize();
    cudaEventRecord(start, 0);
    for (int i = 0; i < 10; i++)
    {
        init1 << <blocksPerGrid, threadsPerBlock >> > (sizes[0], dev_diagCopy, dev_reciprocalDiagPtr, dev_dependenciesCopy);
        init2 << <blocksPerGrid, threadsPerBlock >> > (dev_reciprocalDiagPtr, dev_upper, dev_lowerAddr, dev_upperAddr, dev_dependenciesCopy, dev_dependencies, sizes[1]);
        init3 << <blocksPerGrid, threadsPerBlock >> > (sizes[0], dev_reciprocalDiagPtr);
    }
    cudaEventRecord(stop, 0);

    // Synchronize GPU and CPU
    cudaEventSynchronize(stop);

    // Get time in millis
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Tiempo init %f s\n", (elapsedTime / 1000) / 10);

    //Preconditioner
    cudaDeviceSynchronize();
    cudaEventRecord(start, 0);
    for (int i = 0; i < 10; i++)
    {
        precondition1 << <blocksPerGrid, threadsPerBlock >> > (dev_diag, sizes[0], dev_initialResidual, dev_reciprocalDiagPtr, dev_dependenciesCopy, dev_dependenciesReverseCopy);
        precondition2 << <blocksPerGrid, threadsPerBlock >> > (dev_reciprocalDiagPtr, dev_diag, dev_upper, dev_lowerAddr, dev_upperAddr, dev_dependenciesCopy, dev_dependencies, sizes[1]);
        precondition3 << <blocksPerGrid, threadsPerBlock >> > (dev_reciprocalDiagPtr, dev_diag, dev_upper, dev_lowerAddr, dev_upperAddr, dev_dependenciesReverseCopy, dev_dependenciesReverse, sizes[1]);
    }
    cudaEventRecord(stop, 0);

    // Synchronize GPU and CPU
    cudaEventSynchronize(stop);

    // Get time in millis
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Tiempo precondition %f s\n", (elapsedTime / 1000) / 10);

    // Copy back results
    cudaMemcpy(diagCopy, dev_diag, sizes[0] * sizeof(double), cudaMemcpyDeviceToHost);

    // Check results
    /*for (int j = 0; j < sizes[0]; j++) {
        //cout << j << endl;
        if (fabs(diag[j] - diagCopy[j]) > 1e-12) {
            cout << "index: " << j << " " << diag[j] << " " << diagCopy[j] << endl;
            break;
        }
    }*/

    // Free streams and events
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // Free GPU mem
    cudaFree(dev_diag);
    cudaFree(dev_upper);
    cudaFree(dev_diagCopy);
    cudaFree(dev_upperAddr);
    cudaFree(dev_lowerAddr);
    cudaFree(dev_initialResidual);
    cudaFree(dev_reciprocalDiagPtr);
    cudaFree(dev_dependencies);
    cudaFree(dev_dependenciesCopy);

    //Destroy matrix
    destroy(sizes, diag, lower, upper, lowerAddr, upperAddr, ownerStartAddr, initialResidual, psi, dependencies, dependenciesReverse);

    return 0;
}