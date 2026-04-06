#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <ctime>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
using namespace std;

// Version Order Dependencies
__global__ void init1
(
    const int nCells,
    const double* diagPtr,
    double* reciprocalDiagPtr
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nCells) {
        reciprocalDiagPtr[cell] = diagPtr[cell];
    }
}

__global__ void init2
(
    const int start,
    const int end,
    const int* times,
    const int* upperAddr,
    const int* lowerAddr,
    const double* upper,
    const double* diagPtr,
    double* reciprocalDiagPtr
)
{
    int time = blockDim.x * blockIdx.x + threadIdx.x + start;

    if (time < end) {
        for (int face = times[time]; face < times[time + 1]; time++) {
            reciprocalDiagPtr[upperAddr[time]] -= upper[face] * upper[face] / reciprocalDiagPtr[lowerAddr[face]];
        }
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
    const int* times,
    const int* writeMap,
    const int* readMap,
    const double* upper,
    const double* reciprocalDiagPtr
)
{
    int time = blockDim.x * blockIdx.x + threadIdx.x + start;

    if (time < end) {
        for (int face = times[time]; face < times[time + 1]; face++) {
            diag[writeMap[time]] -= reciprocalDiagPtr[writeMap[time]] * upper[face] * diag[readMap[face]];
        }
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

void destroy
(
    int* sizes,
    double* diag,
    double* initialResidual,
    int* readMap,
    int* writeMap,
    int* newLowerAddr,
    int* newUpperAddr,
    double* newUpperInit,
    double* newUpperPrecondition,
    int* timesInit,
    int* timesIndicesInit,
    int* timesPrecondition,
    int* timesIndicesPrecondition
)
{
    delete[] sizes;
    delete[] diag;
    delete[] initialResidual;
    delete[] readMap;
    delete[] writeMap;
    delete[] newLowerAddr;
    delete[] newUpperAddr;
    delete[] newUpperInit;
    delete[] newUpperPrecondition;
    delete[] timesInit;
    delete[] timesIndicesInit;
    delete[] timesPrecondition;
    delete[] timesIndicesPrecondition;
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

void orderByWriteZone
(
    const int nFaces,
    const int timesIndicesSize,
    const int* timesIndices,
    int** writeMap,
    int** readMap,
    double** upper
)
{
    int i;
    int j;
    int k;
    int size;
    int sortedIndex;
    int* indices;
    int* sortedWriteMap;
    int* sortedReadMap;
    double* sortedUpper;

    sortedWriteMap = new int[nFaces];
    sortedReadMap = new int[nFaces];
    sortedUpper = new double[nFaces];

    for (i = 0; i < timesIndicesSize; i++)
    {
        size = timesIndices[i + 1] - timesIndices[i];
        indices = new int[size];

        for (k = 0, j = timesIndices[i]; j < timesIndices[i + 1]; k++, j++) indices[k] = j;

        sort(indices, indices + size, [&](int a, int b) {
            return (*writeMap)[a] < (*writeMap)[b];
            });

        for (k = 0, j = timesIndices[i]; j < timesIndices[i + 1]; k++, j++)
        {
            sortedWriteMap[j] = (*writeMap)[indices[k]];
            sortedReadMap[j] = (*readMap)[indices[k]];
            sortedUpper[j] = (*upper)[indices[k]];
        }

        delete[](indices);
    }

    delete[]((*writeMap));
    delete[]((*readMap));
    delete[]((*upper));

    (*writeMap) = sortedWriteMap;
    (*readMap) = sortedReadMap;
    (*upper) = sortedUpper;
}

void groupBySameWriteMapIndex
(
    const int nFaces,
    int* timesIndices,
    int** times,
    int** writeMap
)
{
    int* fixedTimes;
    int* fixedWriteMap;
    int indexTime = 1;
    int time = 1;

    (*times)[0] = 0;
    timesIndices[0] = 0;

    for (int i = 1; i < nFaces; i++)
    {
        if (timesIndices[indexTime] == i)
        {
            (*times)[time] = i;
            timesIndices[indexTime] = time;
            time++;
            indexTime++;
            (*writeMap)[time] = (*writeMap)[i];
        }
        else if ((*writeMap)[i] != (*writeMap)[i - 1])
        {
            (*times)[time] = i;
            time++;
            (*writeMap)[time] = (*writeMap)[i];
        }
    }
    timesIndices[indexTime] = time;
    (*times)[time] = nFaces;

    fixedWriteMap = new int[time];
    copy((*writeMap), (*writeMap) + time, fixedWriteMap);

    time++;

    fixedTimes = new int[time];
    copy((*times), (*times) + time, fixedTimes);

    delete[]((*times));
    delete[]((*writeMap));
    (*times) = fixedTimes;
    (*writeMap) = fixedWriteMap;
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
    int** writeMap,
    int** readMap,
    int** timesIndicesInit,
    int** timesInit,
    int** timesIndicesPrecondition,
    int** timesPrecondition,
    double** newUpperInit,
    double** newUpperPrecondition,
    int* parallelSectionsInit,
    int* parallelSectionsPrecondition
)
{
    const int nFaces2 = nFaces * 2;
    int* dependencies;
    int* startTime;
    int* startTimeReverse;
    int* countTimesInit;
    int* countTimesPrecondition;
    int* lastWrite;

    // Init variables    
    dependencies = new int[nFaces];
    startTime = new int[nFaces];
    startTimeReverse = new int[nFaces];
    lastWrite = new int[nCells];

    (*timesInit) = new int[nFaces];
    (*timesPrecondition) = new int[nFaces];
    (*parallelSectionsInit) = 0;
    (*parallelSectionsPrecondition) = 0;

    (*newUpperInit) = new double[nFaces];
    (*newUpperPrecondition) = new double[nFaces2];
    (*newUpperAddr) = new int[nFaces];
    (*newLowerAddr) = new int[nFaces];
    (*readMap) = new int[nFaces2];
    (*writeMap) = new int[nFaces2];

    memset(startTime, 0, nFaces * sizeof(int));
    memset(startTimeReverse, 0, nFaces * sizeof(int));
    memset(lastWrite, -1, nCells * sizeof(int));

    // Define at which time each face can be proccess read
    for (int i = 0; i < nFaces; i++)
    {
        if (startTime[ownerStartAddr[upperAddr[i]]] < 1 + startTime[i])
        {
            for (int j = ownerStartAddr[upperAddr[i]]; j < ownerStartAddr[upperAddr[i] + 1]; j++)
                startTime[j] = 1 + startTime[i];
        }

        if (lastWrite[upperAddr[i]] < startTime[i])
            lastWrite[upperAddr[i]] = startTime[i];


        if (startTime[i] > (*parallelSectionsInit))
            (*parallelSectionsInit) = startTime[i];
    }
    (*parallelSectionsInit)++;

    // Redefine at which time it can be process with last cell write time 
    for (int i = 0; i < nFaces; i++)
        startTime[i] = lastWrite[upperAddr[i]];

    // Define first time to reverse part with last written times of non-reverse side
    for (int i = 0; i < nFaces; i++)
        startTimeReverse[i] = startTime[i] + 1;

    // Define start time at reverse side
    for (int i = 0; i < nFaces; i++)
    {
        if (upperAddr[i] + 1 < nCells && ownerStartAddr[upperAddr[i] + 1] - ownerStartAddr[upperAddr[i]] > 0)
            dependencies[i] = upperAddr[i];
        else dependencies[i] = -1;
    }

    for (int i = nFaces - 1; i >= 0; i--)
    {
        if (dependencies[i] != -1)
            startTimeReverse[i] = 1 + lastWrite[dependencies[i]];
        //else startTimeReverse[i] = 0;

        if (lastWrite[lowerAddr[i]] < startTimeReverse[i])
            lastWrite[lowerAddr[i]] = startTimeReverse[i];

        if (startTimeReverse[i] > (*parallelSectionsPrecondition))
            (*parallelSectionsPrecondition) = startTimeReverse[i];
    }
    (*parallelSectionsPrecondition)++;

    // Redefine at which time it can be process with last cell write time (reverse side)
    for (int i = 0; i < nFaces; i++)
        startTimeReverse[i] = lastWrite[lowerAddr[i]];


    // Init Variables with same length as number of parallel sections
    (*timesIndicesInit) = new int[(*parallelSectionsInit) + 1];
    (*timesIndicesPrecondition) = new int[(*parallelSectionsPrecondition) + 1];

    countTimesInit = new int[(*parallelSectionsInit)];
    countTimesPrecondition = new int[(*parallelSectionsPrecondition)];
    memset(countTimesInit, 0, (*parallelSectionsInit) * sizeof(int));
    memset(countTimesPrecondition, 0, (*parallelSectionsPrecondition) * sizeof(int));

    // Define execution times
    for (int i = 0; i < nFaces; i++)
    {
        countTimesInit[startTime[i]]++;
        countTimesPrecondition[startTime[i]]++;
        countTimesPrecondition[startTimeReverse[i]]++;
    }

    // Define parallel sections for Init
    (*timesIndicesInit)[0] = 0;
    for (int i = 0; i < (*parallelSectionsInit); i++)
        (*timesIndicesInit)[i + 1] = (*timesIndicesInit)[i] + countTimesInit[i];


    // Define parallel sections for Precondition
    (*timesIndicesPrecondition)[0] = 0;
    for (int i = 0; i < (*parallelSectionsPrecondition); i++)
        (*timesIndicesPrecondition)[i + 1] = (*timesIndicesPrecondition)[i] + countTimesPrecondition[i];

    cout << "Size per time init:" << endl;
    for (int i = 0; i < (*parallelSectionsInit); i++)
        cout << (*timesIndicesInit)[i + 1] - (*timesIndicesInit)[i] << " ";
    cout << endl;

    cout << "Size per time precondition:" << endl;
    for (int i = 0; i < (*parallelSectionsPrecondition); i++)
        cout << (*timesIndicesPrecondition)[i + 1] - (*timesIndicesPrecondition)[i] << " ";
    cout << endl;

    // Define new upper and lower sections for init
    for (int i = 0; i < nFaces; i++)
    {
        countTimesInit[startTime[i]]--;
        int time = (*timesIndicesInit)[startTime[i]] + countTimesInit[startTime[i]];
        (*newUpperInit)[time] = upper[i];
        (*newLowerAddr)[time] = lowerAddr[i];
        (*newUpperAddr)[time] = upperAddr[i];
    }

    // Define readMap and writeMap for precondition
    for (int i = 0; i < nFaces; i++)
    {
        // Normal direction
        countTimesPrecondition[startTime[i]]--;
        int time = (*timesIndicesPrecondition)[startTime[i]] + countTimesPrecondition[startTime[i]];
        (*newUpperPrecondition)[time] = upper[i];
        (*readMap)[time] = lowerAddr[i];
        (*writeMap)[time] = upperAddr[i];

        // Reverse direction
        countTimesPrecondition[startTimeReverse[i]]--;
        time = (*timesIndicesPrecondition)[startTimeReverse[i]] + countTimesPrecondition[startTimeReverse[i]];
        (*newUpperPrecondition)[time] = upper[i];
        (*readMap)[time] = upperAddr[i];
        (*writeMap)[time] = lowerAddr[i];
    }

    // Reorder by writeMap (precondition) and by newUpperAddr (init) to speedup with cache
    orderByWriteZone(nFaces, (*parallelSectionsInit), (*timesIndicesInit), newUpperAddr, newLowerAddr, newUpperInit);
    orderByWriteZone(nFaces2, (*parallelSectionsPrecondition), (*timesIndicesPrecondition), writeMap, readMap, newUpperPrecondition);

    // Group by same writeMap (precondition) and by same newUpperAddr index to avoid atomic directive use
    groupBySameWriteMapIndex(nFaces, (*timesIndicesInit), timesInit, newUpperAddr);
    groupBySameWriteMapIndex(nFaces2, (*timesIndicesPrecondition), timesPrecondition, writeMap);

    delete[](startTime);
    delete[](lastWrite);
    delete[](dependencies);
    delete[](countTimesInit);
    delete[](startTimeReverse);
    delete[](countTimesPrecondition);
}

int main() {
    // CPU variables
    int* sizes;
    double* psi;
    double* diag;
    double* lower;
    double* upper;
    int* ownerStartAddr;
    double* initialResidual;

    int* lowerAddr;
    int* upperAddr;

    int* readMap;
    int* writeMap;
    int* newLowerAddr;
    int* newUpperAddr;
    double* newUpperInit;
    double* newUpperPrecondition;

    int* timesInit;
    int* timesIndicesInit;
    int* timesPrecondition;
    int* timesIndicesPrecondition;

    int parallelSectionsInit;
    int parallelSectionsPrecondition;
    // GPU variables
    int* dev_upperAddr;
    int* dev_lowerAddr;
    int* dev_writeMap;
    int* dev_readMap;
    int* dev_timesInit;
    int* dev_timesPrecondition;
    double* dev_diag;
    double* dev_diagCopy;
    double* dev_upperInit;
    double* dev_upperPrecondition;
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
    orderDependencies(sizes[0], sizes[1], upper, upperAddr, lowerAddr, ownerStartAddr, &newUpperAddr, 
                      &newLowerAddr, &writeMap, &readMap, &timesIndicesInit,  &timesInit, &timesIndicesPrecondition, 
                      &timesPrecondition, &newUpperInit, &newUpperPrecondition, &parallelSectionsInit, &parallelSectionsPrecondition);

    delete[] psi;
    delete[] lower;
    delete[] upper;
    delete[] ownerStartAddr;
    delete[] lowerAddr;
    delete[] upperAddr;
    
    // Create events and streams
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Reserve GPU mem
    cudaMalloc((void**)&dev_upperAddr, timesIndicesInit[parallelSectionsInit] * sizeof(int));
    cudaMalloc((void**)&dev_lowerAddr, sizes[3] * sizeof(int));
    cudaMalloc((void**)&dev_writeMap, sizes[3] * 2 * sizeof(int));
    cudaMalloc((void**)&dev_readMap, sizes[3] * 2 * sizeof(int));
    cudaMalloc((void**)&dev_upperInit, sizes[3] * sizeof(double));
    cudaMalloc((void**)&dev_upperPrecondition, sizes[3] * sizeof(double));
    cudaMalloc((void**)&dev_initialResidual, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_diag, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_diagCopy, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_reciprocalDiagPtr, sizes[0] * sizeof(double));
    cudaMalloc((void**)&dev_timesInit, (timesIndicesInit[parallelSectionsInit] + 1) * sizeof(int));
    cudaMalloc((void**)&dev_timesPrecondition, (timesIndicesPrecondition[parallelSectionsPrecondition] + 1) * sizeof(int));

    // Send data to GPU
    cudaMemcpy(dev_upperAddr, newUpperAddr, timesIndicesInit[parallelSectionsInit] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lowerAddr, newLowerAddr, sizes[3] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_writeMap, writeMap, timesIndicesPrecondition[parallelSectionsPrecondition] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_readMap, readMap, sizes[3] * 2 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_timesInit, timesInit, (timesIndicesInit[parallelSectionsInit] + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_timesPrecondition, timesPrecondition, (timesIndicesPrecondition[parallelSectionsPrecondition] + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_upperInit, newUpperInit, sizes[3] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_upperPrecondition, newUpperPrecondition, sizes[3] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_initialResidual, initialResidual, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_diag, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_diagCopy, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);

    // Secuencial option
    //diagCopy = new double[sizes[0]];
    /*init(sizes[0], sizes[3], upper, upperAddr, lowerAddr, diag, &reciprocalDiagPtr);
    precondition(diag, sizes[0], sizes[3], upper, upperAddr, lowerAddr, initialResidual, reciprocalDiagPtr);*/

    // Define number of blocks and threads
    int threadsPerBlock = 64;
    int blocksPerGrid = (sizes[0] + threadsPerBlock - 1) / threadsPerBlock;
    int blocksPerKernel;

    // Init
    cudaEventRecord(start, 0);
    for (int i = 0; i < 1000; i++)
    {
        init1 << <blocksPerGrid, threadsPerBlock >> > (sizes[0], dev_diagCopy, dev_reciprocalDiagPtr);
        for (int j = 0; j < parallelSectionsInit; j++) {
            blocksPerKernel = (timesIndicesInit[j + 1] - timesIndicesInit[j] + threadsPerBlock - 1) / threadsPerBlock;
            init2 << <blocksPerKernel, threadsPerBlock >> > (timesIndicesInit[j], timesIndicesInit[j + 1], dev_timesInit,
                                                             dev_upperAddr, dev_lowerAddr, dev_upperInit, dev_diag, dev_reciprocalDiagPtr);
        }
        init3 << <blocksPerGrid, threadsPerBlock >> > (sizes[0], dev_reciprocalDiagPtr);
    }
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);

    // Synchronize GPU and CPU
    cudaEventSynchronize(stop);

    // Get time in seconds
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time init %f s\n", (elapsedTime / 1000) / 1000);

    // Preconditioner
    /*cudaEventRecord(start, 0);
    for (int i = 0; i < 10000; i++)
    {
        precondition1 << <blocksPerGrid, threadsPerBlock >> > (dev_diag, sizes[0], dev_initialResidual, dev_reciprocalDiagPtr);
        for (int j = 0; j < parallelSectionsPrecondition; j++) {
            blocksPerKernel = (timesIndicesPrecondition[j + 1] - timesIndicesPrecondition[j] + threadsPerBlock - 1) / threadsPerBlock;
            precondition2 << <blocksPerKernel, threadsPerBlock >> > (dev_diag, timesIndicesPrecondition[j], timesIndicesPrecondition[j + 1],
                                                                     dev_timesInit, dev_writeMap, dev_readMap, dev_upperPrecondition, dev_reciprocalDiagPtr);
        }
    }
    cudaEventRecord(stop, 0);

    // Synchronize GPU and CPU
    cudaEventSynchronize(stop);

    // Get time in seconds
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time precondition %f s\n", (elapsedTime / 1000) / 10000);*/

    // Copy back results
    //cudaMemcpy(diagCopy, dev_diag, sizes[0] * sizeof(double), cudaMemcpyDeviceToHost);

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
    cudaFree(dev_readMap);
    cudaFree(dev_writeMap);
    cudaFree(dev_diagCopy);
    cudaFree(dev_upperAddr);
    cudaFree(dev_lowerAddr);
    cudaFree(dev_upperInit);
    cudaFree(dev_timesInit);
    cudaFree(dev_initialResidual);
    cudaFree(dev_timesPrecondition);
    cudaFree(dev_upperPrecondition);
    cudaFree(dev_reciprocalDiagPtr);

    //Destroy matrix
    destroy(sizes, diag, initialResidual, readMap, writeMap, newLowerAddr, newUpperAddr, newUpperInit,
            newUpperPrecondition, timesInit, timesIndicesInit, timesPrecondition, timesIndicesPrecondition);

    return 0;
}