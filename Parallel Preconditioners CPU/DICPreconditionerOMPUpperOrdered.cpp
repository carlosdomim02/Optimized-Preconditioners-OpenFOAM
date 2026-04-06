#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <omp.h>
using namespace std;

#define REPS 100000

void init
(
    const int timesIndicesEnd,
    const int nCells,
    const int nFaces,
    const double* __restrict__ upper,
    const int* __restrict__ upperAddr,
    const int* __restrict__ lowerAddr,
    const int* __restrict__ timesIndices,
    const int* __restrict__ times,
    const double* __restrict__ diagPtr,
    double* __restrict__ reciprocalDiagPtr
)
{    
    #pragma omp parallel
    {       
      int end1;
      int end2;
      int cell;
      int face;
      int time = 0;
      int indexTime;
           
      #pragma omp for simd
      for (cell=0; cell<nCells; cell++)
      {
          reciprocalDiagPtr[cell] = diagPtr[cell];
      }
       
      for(indexTime = 0; indexTime < timesIndicesEnd; indexTime++)
      {
          end1 = timesIndices[indexTime + 1];
          #pragma omp for 
          for (time = timesIndices[indexTime]; time < end1; time++)
          {
              cell = upperAddr[time];
              end2 = times[time + 1];
              #pragma omp simd
              for(face = times[time]; face < end2; face++)
              {
                  reciprocalDiagPtr[cell] -= upper[face] * upper[face] / reciprocalDiagPtr[lowerAddr[face]];
              }
          }         
      }
      
      #pragma omp for simd
      for (cell=0; cell<nCells; cell++)
      {
          reciprocalDiagPtr[cell] = 1.0/reciprocalDiagPtr[cell];
      }
    }
}

void precondition
(
    double* __restrict__ diagPtr,
    const int timesIndicesEnd,
    const int nCells,
    const double* __restrict__ upper,
    const int* __restrict__ writeMap,
    const int* __restrict__ readMap,
    const int* __restrict__ times,
    const int* __restrict__ timesIndices,
    const double* __restrict__ initialResidual,
    const double* __restrict__ reciprocalDiagPtr
)
{
    #pragma omp parallel 
    {
      int end1;
      int end2;
      int cell;
      int face;
      int time = 0;
      int indexTime;
      double tempValue;
      
      #pragma omp for simd
      for (cell=0; cell<nCells; cell++)
      {
          diagPtr[cell] = reciprocalDiagPtr[cell] * initialResidual[cell];
      }
    
      for(indexTime = 0; indexTime < timesIndicesEnd; indexTime++) 
      {
          end1 = timesIndices[indexTime + 1];
          #pragma omp for
          for (time = timesIndices[indexTime]; time < end1; time++)
          {
              tempValue = 0.0;
              cell = writeMap[time];
              end2 = times[time + 1];
              #pragma omp simd
              for(face = times[time]; face < end2; face++)
              {
                  tempValue += upper[face]*diagPtr[readMap[face]];
              }
              diagPtr[cell] -= reciprocalDiagPtr[cell]*tempValue;              
          }                    
      }
    }
}

void destroy
(
    int *sizes,
    double *diag,
    double *initialResidual,
    int *readMap,
    int *writeMap,
    int *newLowerAddr,
    int *newUpperAddr,
    double *newUpperInit,
    double *newUpperPrecondition,
    int *timesInit,
    int *timesIndicesInit,
    int *timesPrecondition,
    int *timesIndicesPrecondition,
    int *losortAddr,
    int *losortStartAddr
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
    delete[] losortAddr;
    delete[] losortStartAddr;
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
        
    for(i = 0; i < timesIndicesSize; i++) 
    {        
        size = timesIndices[i + 1] - timesIndices[i];
        indices = new int[size];
        
        for(k = 0, j = timesIndices[i]; j < timesIndices[i + 1]; k++, j++) indices[k] = j;
    
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
        
    for(int i = 1; i < nFaces; i++)
    {
        if(timesIndices[indexTime] == i)
        {
            (*times)[time] = i;
            timesIndices[indexTime] = time;          
            (*writeMap)[time] = (*writeMap)[i]; 
            time++;
            indexTime++;
        }
        else if((*writeMap)[i] != (*writeMap)[i - 1])
        {
            (*times)[time] = i;  
            (*writeMap)[time] = (*writeMap)[i]; 
            time++;
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
    memset(startTimeReverse, 1, nFaces * sizeof(int));
    memset(lastWrite, -1, nCells * sizeof(int));

    // Define at which time each face can be proccess read
    for(int i = 0; i < nFaces; i++)
    {
        if(startTime[ownerStartAddr[upperAddr[i]]] < 1 + startTime[i])
        {
            for(int j = ownerStartAddr[upperAddr[i]]; j < ownerStartAddr[upperAddr[i] + 1]; j++)
                startTime[j] = 1 + startTime[i];
        }
        
        // Group all wrtting index at same parallel section
        if(lastWrite[upperAddr[i]] < startTime[i])
            lastWrite[upperAddr[i]] = startTime[i];
        
        if(startTime[i] > (*parallelSectionsInit))
            (*parallelSectionsInit) = startTime[i];
    }
    (*parallelSectionsInit)++;
    
    // Redefine at which time it can be process with last cell write time 
    for(int i = 0; i < nFaces; i++)
        startTime[i] = lastWrite[upperAddr[i]];
        
    // Define last writing time for reverse with last reading time of normal mode
    for(int i = 0; i < nCells; i++)
        if(ownerStartAddr[i+1] - ownerStartAddr[i] > 0)
            lastWrite[i] = maxTime(startTime, ownerStartAddr[i], ownerStartAddr[i+1]);

    // Define start time at reverse side    
    for(int i = nFaces - 1; i >= 0; i--)
    {        
        startTimeReverse[i] = 1 + lastWrite[upperAddr[i]]; 
        
        if(lastWrite[lowerAddr[i]] < startTimeReverse[i])
            lastWrite[lowerAddr[i]] = startTimeReverse[i];
            
        if(startTimeReverse[i] > (*parallelSectionsPrecondition))
            (*parallelSectionsPrecondition) = startTimeReverse[i];
    }
    (*parallelSectionsPrecondition)++;
    
    // Redefine at which time it can be process with last cell write time (reverse side)
    for(int i = 0; i < nFaces; i++)
        startTimeReverse[i] = lastWrite[lowerAddr[i]];
    
    // Init Variables with same length as number of parallel sections
    (*timesIndicesInit) = new int[(*parallelSectionsInit) + 1];
    (*timesIndicesPrecondition) = new int[(*parallelSectionsPrecondition) + 1];
    
    countTimesInit = new int[(*parallelSectionsInit)];
    countTimesPrecondition = new int[(*parallelSectionsPrecondition)];
    memset(countTimesInit, 0, (*parallelSectionsInit) * sizeof(int));
    memset(countTimesPrecondition, 0, (*parallelSectionsPrecondition) * sizeof(int));
    
    // Define execution times
    for(int i = 0; i < nFaces; i++)
    {
        countTimesInit[startTime[i]]++; 
        countTimesPrecondition[startTime[i]]++; 
        countTimesPrecondition[startTimeReverse[i]]++; 
    } 
    
    // Define parallel sections for Init
    (*timesIndicesInit)[0] = 0;
    for(int i = 0; i < (*parallelSectionsInit); i++)
        (*timesIndicesInit)[i + 1] = (*timesIndicesInit)[i] + countTimesInit[i];
    
    // Define parallel sections for Precondition
    (*timesIndicesPrecondition)[0] = 0;
    for(int i = 0; i < (*parallelSectionsPrecondition); i++)
        (*timesIndicesPrecondition)[i + 1] = (*timesIndicesPrecondition)[i] + countTimesPrecondition[i];
         
    // Define new upper and lower sections for init
    for(int i = 0; i < nFaces; i++)
    {
        countTimesInit[startTime[i]]--;
        int time = (*timesIndicesInit)[startTime[i]] + countTimesInit[startTime[i]];
        (*newUpperInit)[time] = upper[i];
        (*newLowerAddr)[time] = lowerAddr[i];
        (*newUpperAddr)[time] = upperAddr[i];
    }
    
    // Define readMap and writeMap for precondition
    for(int i = 0; i < nFaces; i++)
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

int main(int argc, char *argv[]){
    int *sizes;
    double *psi;  
    double *diag;
    double *diagCopy;
    double *lower;
    double *upper;
    int *ownerStartAddr;
    int *losortStartAddr;
    double *initialResidual;
    double *reciprocalDiagPtr;
    double *reciprocalDiagPtrCopy;
    
    int *lowerAddr;
    int *upperAddr;
    int *losortAddr;
    
    int *readMap;
    int *writeMap;
    int *newLowerAddr;
    int *newUpperAddr;
    double *newUpperInit;
    double *newUpperPrecondition;
    
    int *timesInit;
    int *timesIndicesInit;
    int *timesPrecondition;
    int *timesIndicesPrecondition;
    
    int parallelSectionsInit;
    int parallelSectionsPrecondition;
      
    double t0, t1;
    int numThreads;
    
    if(argc != 2)
    {
        fprintf(stderr, "Use: need a number of threads <int>\n");
        return -1;
    }
    
    numThreads = atoi(argv[1]);
    
    omp_set_num_threads(numThreads);

    //Read matrix values
    binaryRead("./matrixBeforeLosort.bin", &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi, &ownerStartAddr, &losortAddr, &losortStartAddr);

    //Order dependencies 
    orderDependencies(sizes[0], sizes[1], upper, upperAddr, lowerAddr, ownerStartAddr, &newUpperAddr, &newLowerAddr, &writeMap, &readMap, &timesIndicesInit, 
                      &timesInit, &timesIndicesPrecondition, &timesPrecondition, &newUpperInit, &newUpperPrecondition, &parallelSectionsInit, &parallelSectionsPrecondition);
                      
    delete[] psi;
    delete[] lower;
    delete[] ownerStartAddr;
    
    // Copies to check results
    diagCopy = new double[sizes[0]];
    for (int cell=0; cell<sizes[0]; cell++)
    {
        diagCopy[cell] = diag[cell];
    }
    
    reciprocalDiagPtr = new double[sizes[0]];
    reciprocalDiagPtrCopy = new double[sizes[0]];
    
    // Init
    t0 = omp_get_wtime();	
    for(int i = 0; i < REPS; i++){
        init(parallelSectionsInit, sizes[0], sizes[3], newUpperInit, newUpperAddr, newLowerAddr, timesIndicesInit, timesInit, diag, reciprocalDiagPtr);   
    }
    t1 = omp_get_wtime();	
    
    cout << (t1 - t0) / REPS << endl;
    
    // Precondition
    t0 = omp_get_wtime();	
    for(int i = 0; i < REPS; i++){
        precondition(diag, parallelSectionsPrecondition, sizes[0], newUpperPrecondition, writeMap, readMap, 
                           timesPrecondition, timesIndicesPrecondition, initialResidual, reciprocalDiagPtr);
    }
    t1 = omp_get_wtime();	
    
    cout << (t1 - t0) / REPS << endl;
    
    //cout << "Threads: " << omp_get_max_threads() << endl;
    
    //Check results   
    //Init 
    for (int cell=0; cell<sizes[0]; cell++)
    {
        reciprocalDiagPtrCopy[cell] = diagCopy[cell];
    }
    
    for (int face=0; face<sizes[1]; face++)
    {
        reciprocalDiagPtrCopy[upperAddr[face]] -= upper[face] * upper[face] / reciprocalDiagPtrCopy[lowerAddr[face]];
    }
    
    for (int cell=0; cell<sizes[0]; cell++)
    {
        reciprocalDiagPtrCopy[cell] = 1.0/reciprocalDiagPtrCopy[cell];
    }
    
    //Precondition
    int nFacesM1 = sizes[1] - 1;

    for (int cell=0; cell<sizes[0]; cell++)
    {
        diagCopy[cell] = reciprocalDiagPtrCopy[cell] * initialResidual[cell];
    }

    for (int face=0; face<sizes[1]; face++)
    {
        diagCopy[upperAddr[face]] -= reciprocalDiagPtrCopy[upperAddr[face]]*upper[face]*diagCopy[lowerAddr[face]];
    }

    for (int face=nFacesM1; face>=0; face--)
    {
        diagCopy[lowerAddr[face]] -= reciprocalDiagPtrCopy[lowerAddr[face]]*upper[face]*diagCopy[upperAddr[face]];
    }
    
    //Checking
    for(int j = 0; j < sizes[0]; j++)
    {
        
        if(fabs(reciprocalDiagPtr[j] - reciprocalDiagPtrCopy[j]) > 1e-9)
        {
            cout << "Index: " << j << " error in reciprocalDiagPtr == " << reciprocalDiagPtr[j] << " and reciprocalDiagPtrCopy == " << reciprocalDiagPtrCopy[j] << endl;
            break;
        }
        if(fabs(diag[j] - diagCopy[j]) > 1e-9)
        {
            cout << "Index: " << j << " error in diag == " << diag[j] << " and diagCopy == " << diagCopy[j] << endl;
            break;
        }
    }    
    
    delete[] upper;
    delete[] diagCopy;
    delete[] lowerAddr;
    delete[] upperAddr;
    delete[] reciprocalDiagPtr;
    delete[] reciprocalDiagPtrCopy;

    //Destroy matrix
    destroy(sizes, diag, initialResidual, readMap, writeMap, newLowerAddr, newUpperAddr, newUpperInit, 
            newUpperPrecondition, timesInit, timesIndicesInit, timesPrecondition, timesIndicesPrecondition, losortAddr, losortStartAddr);

    return 0;
}