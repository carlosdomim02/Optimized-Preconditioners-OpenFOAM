#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cooperative_groups.h>
using namespace std;

#define REPS 1000

__global__ void init1
(
    const int nCells,
    const double* diagPtr,
    double* reciprocalDiagPtr
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nCells)
    {
        reciprocalDiagPtr[cell] = diagPtr[cell];
    }
}


__global__ void init2
(
    const int sectionStart,
    const int sectionEnd,
    const double* upper,
    const int* upperAddr,
    const int* lowerAddr,
    double* reciprocalDiagPtr
)
{
    int face = sectionStart + blockDim.x * blockIdx.x + threadIdx.x;

    if (face < sectionEnd)
    {
        reciprocalDiagPtr[upperAddr[face]] -= upper[face] * upper[face] / reciprocalDiagPtr[lowerAddr[face]];
    }
}

__global__ void init3
(
    const int nCells,
    double* reciprocalDiagPtr
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nCells)
    {
        reciprocalDiagPtr[cell] = 1.0 / reciprocalDiagPtr[cell];
    }
}

__global__ void precondition1
(

    const int nCells,
    double* diagPtr,
    const double* initialResidual,
    const double* reciprocalDiagPtr
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nCells)
    {
        diagPtr[cell] = reciprocalDiagPtr[cell] * initialResidual[cell];
    }
}

__global__ void precondition2
(
    const int sectionStart,
    const int sectionEnd,
    double* diagPtr,
    const double* upper,
    const int* writeMap,
    const int* readMap,
    const double* reciprocalDiagPtr
)
{
    int face = sectionStart + blockDim.x * blockIdx.x + threadIdx.x;

    if (face < sectionEnd)
    {
        diagPtr[writeMap[face]] -= reciprocalDiagPtr[writeMap[face]] * upper[face] * diagPtr[readMap[face]];
    }
}

void destroy(
    int* sizes,
    double* psi,
    double* diag,
    double* lower,
    double* upper,
    int* ownerStartAddr,
    int* losortStartAddr,
    double* initialResidual,
    double* reciprocalDiagPtr,
    int* lowerAddr,
    int* upperAddr,
    int* losortAddr,
    int* readMap,
    int* writeMap,
    int* newLowerAddr,
    int* newUpperAddr,
    double* newUpperInit,
    double* newUpperPrecondition,
    int* timesIndicesInit,
    int* timesIndicesPrecondition
)
{
    delete[] sizes;
    delete[] psi;
    delete[] diag;
    delete[] lower;
    delete[] upper;
    delete[] ownerStartAddr;
    delete[] losortStartAddr;
    delete[] initialResidual;
    delete[] reciprocalDiagPtr;
    delete[] lowerAddr;
    delete[] upperAddr;
    delete[] losortAddr;
    delete[] readMap;
    delete[] writeMap;
    delete[] newLowerAddr;
    delete[] newUpperAddr;
    delete[] newUpperInit;
    delete[] newUpperPrecondition;
    delete[] timesIndicesInit;
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
    int** writeMap,
    int** readMap,
    int** timesIndicesInit,
    int** timesIndicesPrecondition,
    double** newUpperInit,
    double** newUpperPrecondition,
    int* parallelSectionsInit,
    int* parallelSectionsPrecondition
)
{
    const int nFaces2 = nFaces * 2;
    int* startTime;
    int* startTimeReverse;
    int* countTimesInit;
    int* countTimesPrecondition;
    int* lastWrite;

    // Init variables    
    startTime = new int[nFaces];
    startTimeReverse = new int[nFaces];
    lastWrite = new int[nCells];

    (*parallelSectionsInit) = 0;
    (*parallelSectionsPrecondition) = 0;

    (*newUpperInit) = new double[nFaces];
    (*newUpperPrecondition) = new double[nFaces2];
    (*newUpperAddr) = new int[nFaces];
    (*newLowerAddr) = new int[nFaces];
    (*readMap) = new int[nFaces2];
    (*writeMap) = new int[nFaces2];

    memset(startTime, 0, nFaces * sizeof(int));
    memset(startTimeReverse, 1, nFaces * sizeof(int));
    memset(lastWrite, -1, nCells * sizeof(int));

    // Define at which time each face can be proccess read
    for (int i = 0; i < nFaces; i++)
    {
        // Break up groups of same writting index at the same time
        if (lastWrite[upperAddr[i]] >= startTime[i])
        {
            lastWrite[upperAddr[i]]++;
            startTime[i] = lastWrite[upperAddr[i]];
        }
        else lastWrite[upperAddr[i]] = startTime[i];

        if (startTime[ownerStartAddr[upperAddr[i]]] < 1 + startTime[i])
        {
            for (int j = ownerStartAddr[upperAddr[i]]; j < ownerStartAddr[upperAddr[i] + 1]; j++)
                startTime[j] = 1 + startTime[i];
        }

        if (startTime[i] > (*parallelSectionsInit))
            (*parallelSectionsInit) = startTime[i];
    }
    (*parallelSectionsInit)++;

    // Define last writing time for reverse with last reading time of normal mode
    for (int i = 0; i < nCells; i++)
        if (ownerStartAddr[i + 1] - ownerStartAddr[i] > 0)
            lastWrite[i] = maxTime(startTime, ownerStartAddr[i], ownerStartAddr[i + 1]);
    // else lastWrite[i] = pastLastWrite

    // Define start time at reverse side
    for (int i = nFaces - 1; i >= 0; i--)
    {
        startTimeReverse[i] = 1 + lastWrite[upperAddr[i]];

        if (lastWrite[lowerAddr[i]] >= startTimeReverse[i]) {
            lastWrite[lowerAddr[i]]++;
            startTimeReverse[i] = lastWrite[lowerAddr[i]];
        }
        else lastWrite[lowerAddr[i]] = startTimeReverse[i];

        if (startTimeReverse[i] > (*parallelSectionsPrecondition))
            (*parallelSectionsPrecondition) = startTimeReverse[i];
    }
    (*parallelSectionsPrecondition)++;

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

    // Print counts
    /*cout << "Size per time init:" << endl;
    int countInit = 0, countPre = 0, countRev = 0;
    for (int i = 0; i < (*parallelSectionsInit); i++) {
        cout << (*timesIndicesInit)[i + 1] - (*timesIndicesInit)[i] << " + ";
        countInit += countTimesInit[i];
    }
    cout << endl;

    cout << "Size per time precondition:" << endl;
    for (int i = 0; i < (*parallelSectionsPrecondition); i++) {
        cout << (*timesIndicesPrecondition)[i + 1] - (*timesIndicesPrecondition)[i] << " + ";
    }
    cout << endl;*/

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

    delete[](startTime);
    delete[](lastWrite);
    delete[](countTimesInit);
    delete[](startTimeReverse);
    delete[](countTimesPrecondition);
}

int main(int argc, char* argv[]) {
    // CPU variables
    int* sizes;
    double* psi;
    double* diag;
    double* diagCopy;
    double* lower;
    double* upper;
    int* ownerStartAddr;
    int* losortStartAddr;
    double* initialResidual;
    double* reciprocalDiagPtr;
    double* reciprocalDiagPtrCopy;
    int* lowerAddr;
    int* upperAddr;
    int* losortAddr;
    int* readMap;
    int* writeMap;
    int* newLowerAddr;
    int* newUpperAddr;
    double* newUpperInit;
    double* newUpperPrecondition;
    int* timesIndicesInit;
    int* timesIndicesPrecondition;
    int parallelSectionsInit;
    int parallelSectionsPrecondition;
    // GPU variables
    int* dev_upperAddr; //sizes[1]
    int* dev_lowerAddr; //sizes[1]
    int* dev_writeMap; //sizes[1] * 2
    int* dev_readMap; //sizes[1] * 2
    double* dev_diag; //sizes[0]
    double* dev_upperInit; //sizes[1]
    double* dev_upperPrecondition; //sizes[1] * 2
    double* dev_initialResidual; //sizes[0]
    double* dev_reciprocalDiagPtr; //sizes[0]
    // CUDA time measuring variables
    cudaEvent_t start;
    cudaEvent_t stop;
    float elapsedTime;

    cudaSetDevice(0);

    // Read matrix values
    binaryRead("./matrixBeforeLong.bin", &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi, &ownerStartAddr, &losortAddr, &losortStartAddr);

    // Order dependencies 
    orderDependencies(sizes[0], sizes[1], upper, upperAddr, lowerAddr, ownerStartAddr, &newUpperAddr, &newLowerAddr, &writeMap, &readMap, &timesIndicesInit,
        &timesIndicesPrecondition, &newUpperInit, &newUpperPrecondition, &parallelSectionsInit, &parallelSectionsPrecondition);

    // Create time events
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Reserve GPU mem
    cudaMalloc(&dev_upperAddr, sizes[1] * sizeof(int));
    cudaMalloc(&dev_lowerAddr, sizes[1] * sizeof(int));
    cudaMalloc(&dev_writeMap, sizes[1] * 2 * sizeof(int));
    cudaMalloc(&dev_readMap, sizes[1] * 2 * sizeof(int));
    cudaMalloc(&dev_diag, sizes[0] * sizeof(double));
    cudaMalloc(&dev_upperInit, sizes[1] * sizeof(double));
    cudaMalloc(&dev_upperPrecondition, sizes[1] * 2 * sizeof(double));
    cudaMalloc(&dev_initialResidual, sizes[0] * sizeof(double));
    cudaMalloc(&dev_reciprocalDiagPtr, sizes[0] * sizeof(double));

    // Send data to GPU
    cudaMemcpy(dev_upperAddr, newUpperAddr, sizes[1] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lowerAddr, newLowerAddr, sizes[1] * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_writeMap, writeMap, sizes[1] * 2 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_readMap, readMap, sizes[1] * 2 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_diag, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_upperInit, newUpperInit, sizes[1] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_upperPrecondition, newUpperPrecondition, sizes[1] * 2 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_initialResidual, initialResidual, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);

    // Checking vectors config
    diagCopy = new double[sizes[0]];

    reciprocalDiagPtr = new double[sizes[0]];
    reciprocalDiagPtrCopy = new double[sizes[0]];

    for (int k = 8; k <= 1024; k *= 2)
    {
        int startSection;
        int endSection;
        int threadsPerBlock = k;
        int blocksPerGridSection;
        int blocksPerGridCell = (sizes[0] + threadsPerBlock - 1) / threadsPerBlock;

        cout << "DIC Preconditioner CUDA (" << k << " Num Threads):" << endl;

        for (int j = 0; j < 2; j++)
        {
            cudaMemcpy(dev_diag, diag, sizes[0] * sizeof(double), cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();

            // Time init
            cudaEventRecord(start, 0);
            for (int i = 0; i < REPS; i++) {
                init1 << <blocksPerGridCell, threadsPerBlock >> > (sizes[0], dev_diag, dev_reciprocalDiagPtr);
                cudaDeviceSynchronize();

                startSection = 0;
                endSection = timesIndicesInit[1];
                blocksPerGridSection = (endSection - startSection + threadsPerBlock - 1) / threadsPerBlock;
                for (int i = 1; i <= parallelSectionsInit; i++)
                {
                    init2 << <blocksPerGridSection, threadsPerBlock >> >
                        (startSection, endSection, dev_upperInit, dev_upperAddr, dev_lowerAddr, dev_reciprocalDiagPtr);
                    if (i < parallelSectionsInit)
                    {
                        startSection = timesIndicesInit[i];
                        endSection = timesIndicesInit[i + 1];
                        blocksPerGridSection = (endSection - startSection + threadsPerBlock - 1) / threadsPerBlock;
                    }
                    cudaDeviceSynchronize();
                }

                init3 << <blocksPerGridCell, threadsPerBlock >> > (sizes[0], dev_reciprocalDiagPtr);
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
                precondition1 << <blocksPerGridCell, threadsPerBlock >> > (sizes[0], dev_diag, dev_initialResidual, dev_reciprocalDiagPtr);
                cudaDeviceSynchronize();

                startSection = 0;
                endSection = timesIndicesPrecondition[1];
                blocksPerGridSection = (endSection - startSection + threadsPerBlock - 1) / threadsPerBlock;
                for (int i = 1; i <= parallelSectionsPrecondition; i++)
                {
                    precondition2 << <blocksPerGridSection, threadsPerBlock >> >
                        (startSection, endSection, dev_diag, dev_upperPrecondition, dev_writeMap, dev_readMap, dev_reciprocalDiagPtr);
                    if (i < parallelSectionsPrecondition)
                    {
                        startSection = timesIndicesPrecondition[i];
                        endSection = timesIndicesPrecondition[i + 1];
                        blocksPerGridSection = (endSection - startSection + threadsPerBlock - 1) / threadsPerBlock;
                    }
                    cudaDeviceSynchronize();
                }
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
        for (int i = 0; i < 1000; i++)
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
        cout << (double(t1 - t0) / CLOCKS_PER_SEC) / 1000 << endl;

        // Precondition
        t0 = clock();
        int nFacesM1 = sizes[1] - 1;
        for (int i = 0; i < 1000; i++)
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
        cout << (double(t1 - t0) / CLOCKS_PER_SEC) / 1000 << endl;
    }

    // Return to CPU result data to compare
    cudaMemcpy(diag, dev_diag, sizes[0] * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(reciprocalDiagPtr, dev_reciprocalDiagPtr, sizes[0] * sizeof(double), cudaMemcpyDeviceToHost);

    // Free GPU mem
    cudaFree(dev_upperAddr);
    cudaFree(dev_lowerAddr);
    cudaFree(dev_writeMap);
    cudaFree(dev_readMap);
    cudaFree(dev_diag);
    cudaFree(dev_upperInit);
    cudaFree(dev_upperPrecondition);
    cudaFree(dev_initialResidual); 
    cudaFree(dev_reciprocalDiagPtr);

    // Checking
    for (int j = 0; j < sizes[0]; j++)
    {

        if (fabs(reciprocalDiagPtr[j] - reciprocalDiagPtrCopy[j]) > 1e-12)
        {
            cout << "Index: " << j << " error in reciprocalDiagPtr == " << reciprocalDiagPtr[j] << " and reciprocalDiagPtrCopy == " << reciprocalDiagPtrCopy[j] << endl;
            break;
        }
        if (fabs(diag[j] - diagCopy[j]) > 1e-12)
        {
            cout << "Index: " << j << " error in diag == " << diag[j] << " and diagCopy == " << diagCopy[j] << endl;
            break;
        }
    }

    delete[](diagCopy);
    delete[](reciprocalDiagPtrCopy);

    //Destroy matrix
    destroy(sizes, psi, diag, lower, upper, ownerStartAddr, losortStartAddr, initialResidual,
        reciprocalDiagPtr, lowerAddr, upperAddr, losortAddr, readMap, writeMap, newLowerAddr,
        newUpperAddr, newUpperInit, newUpperPrecondition, timesIndicesInit, timesIndicesPrecondition);

    return 0;
}