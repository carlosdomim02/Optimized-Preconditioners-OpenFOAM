#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cooperative_groups.h>
using namespace std;

#define REPS 100000000

// Version Order Dependencies
__global__ void init1
(
    const int nCells,
    const double* diagPtr,
    double* reciprocalDiagPtr
)
{
    for (int k = 0; k < REPS; k++) {
        int cell = blockDim.x * blockIdx.x + threadIdx.x;

        if (cell < nCells) {
            reciprocalDiagPtr[cell] = diagPtr[cell];
        }
    }
}

/*__global__ void init2
(
    const int start,
    const int end,
    const int* upperAddr,
    const int* lowerAddr,
    const double* upper,
    const double* diagPtr,
    double* reciprocalDiagPtr
)
{
    int face = blockDim.x * blockIdx.x + threadIdx.x + start;

    if (face < end) {
        reciprocalDiagPtr[upperAddr[face]] -= upper[face] * upper[face] / reciprocalDiagPtr[lowerAddr[face]];
    }
}*/

__global__ void init2
(
    const int offsetsSize,
    const int* offsets,
    const int* upperAddr,
    const int* lowerAddr,
    const double* upper,
    const double* diagPtr,
    double* reciprocalDiagPtr
)
{
    for (int k = 0; k < REPS; k++) {
        int i = 0;
        int face;

        while (i < offsetsSize) {
            face = blockDim.x * blockIdx.x + threadIdx.x + offsets[i];
            if (face < offsets[i + 1]) {
                reciprocalDiagPtr[upperAddr[face]] -= upper[face] * upper[face] / reciprocalDiagPtr[lowerAddr[face]];
            }
            i++;
            __syncthreads();
            //cooperative_groups::this_grid().sync();
        }
    }
}

__global__ void init3
(
    const int nCells,
    double* reciprocalDiagPtr
)
{
    for (int k = 0; k < REPS; k++) {
        int cell = blockDim.x * blockIdx.x + threadIdx.x;

        if (cell < nCells) {
            reciprocalDiagPtr[cell] = 1.0 / reciprocalDiagPtr[cell];
        }
    }
}

__global__ void precondition1
(
    double* diag,
    const int nCells,
    const double* initialResidual,
    const double* reciprocalDiagPtr
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nCells) {
        diag[cell] = reciprocalDiagPtr[cell] * initialResidual[cell];
    }
}

__global__ void precondition2
(
    double* diag,
    const int start,
    const int end,
    const int* upperAddr,
    const int* lowerAddr,
    const double* upper,
    const double* reciprocalDiagPtr
)
{
    int face = blockDim.x * blockIdx.x + threadIdx.x + start;

    if (face < end) {
        diag[upperAddr[face]] -= reciprocalDiagPtr[upperAddr[face]] * upper[face] * diag[lowerAddr[face]];
    }
}

__global__ void precondition3
(
    double* diag,
    const int start,
    const int end,
    const int* upperAddr,
    const int* lowerAddr,
    const double* upper,
    const double* reciprocalDiagPtr
)
{
    int face = blockDim.x * blockIdx.x + threadIdx.x + start;

    if (face < end) {
        diag[lowerAddr[face]] -= reciprocalDiagPtr[lowerAddr[face]] * upper[face] * diag[upperAddr[face]];
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
    double* newUpper,
    double* newUpperReverse,
    int* lowerAddr,
    int* upperAddr,
    int* newLowerAddr,
    int* newUpperAddr,
    int* newLowerAddrReverse,
    int* newUpperAddrReverse,
    int* ownerStartAddr,
    int* timesIndices,
    int* timesIndicesReverse,
    double* initialResidual,
    double* psi)
{
    delete[] sizes;
    delete[] diag;
    delete[] lower;
    delete[] upper;
    delete[] newUpper;
    delete[] newUpperReverse;
    delete[] lowerAddr;
    delete[] upperAddr;
    delete[] newLowerAddr;
    delete[] newUpperAddr;
    delete[] newLowerAddrReverse;
    delete[] newUpperAddrReverse;
    delete[] ownerStartAddr;
    delete[] timesIndices;
    delete[] timesIndicesReverse;
    delete[] initialResidual;
    delete[] psi;
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

int maxTime
(
    const int* startTime,
    const int start,
    const int end
)
{
    int max = startTime[start];
    for (int i = start + 1; i < end; i++)
        if (startTime[i] > max)
            max = startTime[i];
    return max;
}

void orderDependencies
(
    const int nCells,
    const int nFaces,
    const double* upper,
    const int* upperAddr,
    const int* lowerAddr,
    const int* ownerStartAddr,
    int** newUpperAddr,
    int** newLowerAddr,
    int** timesIndices,
    double** newUpper,
    int* maxTimes, 
    int* timesIndicesSize
)
{
    int* cellWriteTime;
    int* dependencies;
    int* startTime;
    int* countTimes;

    // Init variables
    (*newUpper) = new double[nFaces];
    startTime = new int[nFaces];
    (*timesIndices) = new int[nCells];
    dependencies = new int[nFaces];
    *newLowerAddr = new int[nFaces];
    *newUpperAddr = new int[nFaces];
    cellWriteTime = new int[nCells];
    countTimes = new int[nCells];
    memset(cellWriteTime, -1, nCells * sizeof(int));
    memset(dependencies, -1, nFaces * sizeof(int));
    memset(startTime, 0, nFaces * sizeof(int));
    memset(countTimes, 0, nCells * sizeof(int));

    for (int i = 0; i < nFaces; i++) {
        // Upgrade dependency in all faces with same owner
        if (cellWriteTime[upperAddr[i]] >= startTime[i]) {
            cellWriteTime[upperAddr[i]]++;
            startTime[i] = cellWriteTime[upperAddr[i]];
        }
        else cellWriteTime[upperAddr[i]] = startTime[i];

        if (startTime[ownerStartAddr[upperAddr[i]]] < 1 + startTime[i])
            for (int j = ownerStartAddr[upperAddr[i]]; j < ownerStartAddr[upperAddr[i] + 1]; j++)
                startTime[j] = 1 + startTime[i];
    }

    for (int i = 0; i < nFaces; i++) {
        countTimes[startTime[i]]++;
    }

    (*timesIndices)[0] = 0;
    for (int i = 0; i < nCells; i++) {
        if (countTimes[i] == 0) {
            (*timesIndices)[i + 1] = 0;
            break;
        }
        (*timesIndices)[i + 1] = (*timesIndices)[i] + countTimes[i];
    }

    // Search max countTimes and timesIndicesSize
    (*maxTimes) = 0;
    (*timesIndicesSize) = 0;
    for (int i = 0; i < nCells; i++) {
        if (countTimes[i] == 0) {
            break;
        }
        if (countTimes[i] > (*maxTimes)) {
            (*maxTimes) = countTimes[i];
        }
        (*timesIndicesSize)++;
    }

    // Define new lower and upper Addrs
    for (int i = 0; i < nFaces; i++) {
        countTimes[startTime[i]]--;
        int time = (*timesIndices)[startTime[i]] + countTimes[startTime[i]];
        (*newUpper)[time] = upper[i];
        (*newLowerAddr)[time] = lowerAddr[i];
        (*newUpperAddr)[time] = upperAddr[i];
    }

    delete[](startTime);
    delete[](countTimes);
    delete[](dependencies);
    delete[](cellWriteTime);
}

void orderDependenciesReverse
(
    const int nCells,
    const int nFaces,
    const double* upper,
    const int* upperAddr,
    const int* lowerAddr,
    const int* ownerStartAddr,
    int** newUpperAddr,
    int** newLowerAddr,
    int** timesIndices,
    double** newUpper
)
{
    int* cellWriteTime;
    int* dependencies;
    int* startTime;
    int* countTimes;

    // Init variables
    (*newUpper) = new double[nFaces];
    startTime = new int[nFaces];
    (*timesIndices) = new int[nCells];
    dependencies = new int[nFaces];
    (*newLowerAddr) = new int[nFaces];
    (*newUpperAddr) = new int[nFaces];
    cellWriteTime = new int[nCells];
    countTimes = new int[nCells];
    memset(cellWriteTime, -1, nCells * sizeof(int));
    memset(countTimes, 0, nCells * sizeof(int));

    // Searching for dependencies
    for (int i = 0; i < nFaces; i++) {
        if (upperAddr[i] + 1 < nCells && ownerStartAddr[upperAddr[i] + 1] - ownerStartAddr[upperAddr[i]] > 0)
            dependencies[i] = upperAddr[i];
        else dependencies[i] = -1;
    }

    // Define execution times
    countTimes[0]++;
    startTime[nFaces - 1] = 0;
    for (int i = nFaces - 2; i >= 0; i--) {
        if (dependencies[i] != -1)
            startTime[i] = 1 + maxTime(startTime, ownerStartAddr[dependencies[i]], ownerStartAddr[dependencies[i] + 1]);
        else startTime[i] = 0;

        if (cellWriteTime[lowerAddr[i]] >= startTime[i]) {
            cellWriteTime[lowerAddr[i]]++;
            startTime[i] = cellWriteTime[lowerAddr[i]];
        }
        else cellWriteTime[lowerAddr[i]] = startTime[i];

        countTimes[startTime[i]]++;
    }

    (*timesIndices)[0] = 0;
    for (int i = 0; i < nCells; i++) {
        if (countTimes[i] == 0) {
            (*timesIndices)[i + 1] = 0;
            break;
        }
        (*timesIndices)[i + 1] = (*timesIndices)[i] + countTimes[i];
    }

    // Define new lower and upper Addrs
    for (int i = 0; i < nFaces; i++) {
        countTimes[startTime[i]]--;
        int time = (*timesIndices)[startTime[i]] + countTimes[startTime[i]];
        (*newUpper)[time] = upper[i];
        (*newLowerAddr)[time] = lowerAddr[i];
        (*newUpperAddr)[time] = upperAddr[i];
    }

    delete[](startTime);
    delete[](countTimes);
    delete[](dependencies);
    delete[](cellWriteTime);
}

int main() {
    // CPU variables
    int* sizes;
    double* diag;
    double* diagCopy;
    double* lower;
    double* upper;
    double* newUpper;
    double* newUpperReverse;
    double* reciprocalDiagPtr;
    int* lowerAddr;
    int* upperAddr;
    int* newLowerAddr;
    int* newUpperAddr;
    int* newLowerAddrReverse;
    int* newUpperAddrReverse;
    int* ownerStartAddr;
    int* timesIndices;
    int* timesIndicesReverse;
    double* initialResidual;
    double* psi;
    int maxTimes;
    int timesIndicesSize;
    // GPU variables
    int* dev_upperAddr;
    int* dev_lowerAddr;
    int* dev_timesIndices;
    int* dev_upperAddrReverse;
    int* dev_lowerAddrReverse;
    double* dev_diag;
    double* dev_diagCopy;
    double* dev_upper;
    double* dev_upperReverse;
    double* dev_initialResidual;
    double* dev_reciprocalDiagPtr;
    // Special cuda variables
    cudaEvent_t start;
    cudaEvent_t stop;
    float elapsedTime = 0.0f;
    // File variables
    ofstream outputFile;
    outputFile.open("output.txt");

    cudaSetDevice(1);

    //Read matrix values
    binaryRead("./matrixBefore.bin", &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi, &ownerStartAddr);
    
    //Ordder Dependencies
    orderDependencies(sizes[0], sizes[3], upper, upperAddr, lowerAddr, ownerStartAddr, &newUpperAddr, &newLowerAddr, &timesIndices, &newUpper, &maxTimes, &timesIndicesSize);
    orderDependenciesReverse(sizes[0], sizes[3], upper, upperAddr, lowerAddr, ownerStartAddr, &newUpperAddrReverse, &newLowerAddrReverse, &timesIndicesReverse, &newUpperReverse);

    cout << maxTimes << " " << timesIndicesSize << endl;

    // Create events and streams
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Reserve GPU mem
    cudaMalloc((void**)&dev_upperAddr, sizes[3] * sizeof(int));
    cudaMalloc((void**)&dev_lowerAddr, sizes[3] * sizeof(int));
    cudaMalloc((void**)&dev_timesIndices, sizes[3] * sizeof(int));
    cudaMalloc((void**)&dev_upperAddrReverse, sizes[3] * sizeof(int));
    cudaMalloc((void**)&dev_lowerAddrReverse, sizes[3] * sizeof(int));
    cudaMalloc((void**)&dev_upper, sizes[3] * sizeof(double));
    cudaMalloc((void**)&dev_upperReverse, sizes[3] * sizeof(double));
    cudaMalloc((void**)&dev_initialResidual, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_diag, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_diagCopy, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_reciprocalDiagPtr, sizes[0] * sizeof(double));

    // Send data to GPU
    cudaMemcpy(dev_upperAddr, newUpperAddr, sizes[3] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lowerAddr, newLowerAddr, sizes[3] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_timesIndices, timesIndices, sizes[0] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_upperAddrReverse, newUpperAddrReverse, sizes[3] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lowerAddrReverse, newLowerAddrReverse, sizes[3] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_upper, newUpper, sizes[3] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_upperReverse, newUpperReverse, sizes[3] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_initialResidual, initialResidual, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_diag, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_diagCopy, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);

    // Secuencial option
    diagCopy = new double[sizes[0]];
    /*init(sizes[0], sizes[3], upper, upperAddr, lowerAddr, diag, &reciprocalDiagPtr);
    precondition(diag, sizes[0], sizes[3], upper, upperAddr, lowerAddr, initialResidual, reciprocalDiagPtr);*/

    // Define number of blocks and threads
    int threadsPerBlock = 16;
    int blocksPerGrid = (sizes[0] + threadsPerBlock - 1) / threadsPerBlock;
    int blocksPerKernel;

    //Preconditioner
    cudaEventRecord(start, 0);
    init1 << <1, sizes[0] >> > (sizes[0], dev_diagCopy, dev_reciprocalDiagPtr);
    cudaDeviceSynchronize();
    //blocksPerKernel = (maxTimes + threadsPerBlock - 1) / threadsPerBlock;
    init2 << <1, maxTimes >> > (timesIndicesSize, dev_timesIndices, dev_upperAddr, dev_lowerAddr, dev_upper, dev_diag, dev_reciprocalDiagPtr);
    cudaDeviceSynchronize();
    init3 << <1, sizes[0] >> > (sizes[0], dev_reciprocalDiagPtr);
    cudaDeviceSynchronize();
        /*outputFile << endl;
        //cudaDeviceSynchronize();
        outputFile << "precondition1 =>\t bloques: " << blocksPerGrid << " hilos: " << threadsPerBlock
            << " trabajo: " << blocksPerGrid * threadsPerBlock << " trabajo real: " << sizes[0] << endl;*/
        //precondition1 << <blocksPerGrid, threadsPerBlock >> > (dev_diag, sizes[0], dev_initialResidual, dev_reciprocalDiagPtr);
        //cudaDeviceSynchronize();
        //outputFile << endl;
        /*for (int j = 0; j < sizes[0]; j++) {
            //cudaDeviceSynchronize();
            if (timesIndices[j + 1] <= 0) break;
            blocksPerKernel = (timesIndices[j + 1] - timesIndices[j] + threadsPerBlock - 1) / threadsPerBlock;
            outputFile << "precondition2 =>\t bloques: " << blocksPerKernel << " hilos: " << threadsPerBlock << " trabajo: "
                       << blocksPerKernel * threadsPerBlock << " trabajo real: "
                       << timesIndices[j + 1] - timesIndices[j] << endl;
            //cudaEventRecord(start, 0);
            precondition2 << <blocksPerKernel, threadsPerBlock >> >
                    (dev_diag, timesIndices[j], timesIndices[j + 1], dev_upperAddr, dev_lowerAddr, dev_upper, dev_reciprocalDiagPtr);
            cudaEventRecord(stop, 0);
            cudaEventSynchronize(stop);
            cudaEventElapsedTime(&elapsedTime, start, stop);
            outputFile << "precondition2 =>\t\t bloques: " << blocksPerKernel << " hilos: " << threadsPerBlock << " trabajo: "
                << blocksPerKernel * threadsPerBlock << " trabajo real: " << timesIndices[j + 1] - timesIndices[j]
                << " tiempo (ms): " << elapsedTime << endl;
        }*/
        //outputFile << endl;
        /*for (int j = 0; j < sizes[0]; j++) {
            //cudaDeviceSynchronize();
            if (timesIndicesReverse[j + 1] <= 0) break;
            blocksPerKernel = (timesIndicesReverse[j + 1] - timesIndicesReverse[j] + threadsPerBlock - 1) / threadsPerBlock;
            outputFile << "precondition3 =>\t bloques: " << blocksPerKernel << " hilos: " << threadsPerBlock << " trabajo: "
                       << blocksPerKernel * threadsPerBlock << " trabajo real: " 
                       << timesIndicesReverse[j + 1] - timesIndicesReverse[j] << endl;
            //cudaEventRecord(start, 0);
            precondition3 << <blocksPerKernel, threadsPerBlock >> >
                (dev_diag, timesIndicesReverse[j], timesIndicesReverse[j + 1], dev_upperAddrReverse, dev_lowerAddrReverse, dev_upperReverse, dev_reciprocalDiagPtr);
            cudaEventRecord(stop, 0);
            cudaEventSynchronize(stop);
            cudaEventElapsedTime(&elapsedTime, start, stop);
            outputFile << "precondition3 =>\t\t bloques: " << blocksPerKernel << " hilos: " << threadsPerBlock << " trabajo: "
                << blocksPerKernel * threadsPerBlock << " trabajo real: " << timesIndices[j + 1] - timesIndices[j]
                << " tiempo (ms): " << elapsedTime << endl;
        }
        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&elapsedTime, start, stop);
        cout << "Time (ms) " << elapsedTime << endl;*/
    cudaEventRecord(stop, 0);

    // Synchronize GPU and CPU
    cudaEventSynchronize(stop);

    // Get time in millis
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Tiempo total %.6ef s\n", (elapsedTime / 1000) / REPS);

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
    cudaFree(dev_upperReverse);
    cudaFree(dev_timesIndices);
    cudaFree(dev_initialResidual);
    cudaFree(dev_upperAddrReverse);
    cudaFree(dev_lowerAddrReverse);
    cudaFree(dev_reciprocalDiagPtr);

    //Destroy matrix
    destroy(sizes, diag, lower, upper, newUpper, newUpperReverse, lowerAddr, upperAddr, newLowerAddr, newUpperAddr, newLowerAddrReverse, newUpperAddrReverse, ownerStartAddr, timesIndices, timesIndicesReverse, initialResidual, psi);

    return 0;
}