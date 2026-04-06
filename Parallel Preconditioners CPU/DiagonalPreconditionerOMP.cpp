#include <immintrin.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstring>
#include <omp.h>
using namespace std;

#define REPS 1000000

void init
(
    const int nCells,
    const double* __restrict__ diagPtr,
    double* __restrict__ reciprocalDiagPtr,
    const double* __restrict__ psi
)
{
    #pragma omp parallel for simd 
    for (int cell=0; cell<nCells; cell++)
    {
        reciprocalDiagPtr[cell] = 1.0/diagPtr[cell];
    }
}

void precondition
(
    double* __restrict__ diag,
    const int nCells,
    const double* __restrict__ initialResidual,
    const double* __restrict__ reciprocalDiagPtr
) 
{
    double t0, t1;
    int sizePerThread;
    int np;
    
    #pragma omp parallel for simd
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
    const double* reciprocalDiagPtr,
    const int* ownerStartAddr
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
    delete[](reciprocalDiagPtr);
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
            FILE *binaryFile = fopen(fileName, "rb");   
            
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

            fread(*diag           , sizeof(double), (*sizes)[0] , binaryFile);
            fread(*lower          , sizeof(double), (*sizes)[1] , binaryFile);
            fread(*lowerAddr      , sizeof(int)   , (*sizes)[2] , binaryFile);
            fread(*upper          , sizeof(double), (*sizes)[3] , binaryFile);
            fread(*upperAddr      , sizeof(int)   , (*sizes)[4] , binaryFile);
            fread(*initialResidual, sizeof(double), (*sizes)[5] , binaryFile);
            fread(*psi            , sizeof(double), (*sizes)[6] , binaryFile);
            fread(*ownerStartAddr , sizeof(int)   , (*sizes)[7] , binaryFile);
          
            fclose(binaryFile);
}

int main(int argc, char *argv[]){
    int *sizes;
    double *diag;
    double *diagCopy;
    double *lower;
    double *upper;
    int *lowerAddr;
    int *upperAddr;
    int *ownerStartAddr;
    double *reciprocalDiagPtrCopy;
    double *reciprocalDiagPtr;
    double *initialResidual;
    double *psi;
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
    binaryRead("../MatricesConXIT300/matrixBefore.bin", &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi, &ownerStartAddr);
    
    //Copy diag for checking phase
    diagCopy = new double[sizes[0]];
    for(int i = 0; i < sizes[0]; i++)
    {
      diagCopy[i] = diag[i];
    }

    reciprocalDiagPtr = new double[sizes[0]];

    //Init
    t0 = omp_get_wtime();	
    for(int k = 0; k < REPS; k++)
    {   
        init(sizes[0], diag, reciprocalDiagPtr, psi);
    }
    t1 = omp_get_wtime();
    cout << "Time init: " << (t1 - t0) / REPS << endl;
    
    //Preconditioner
    t0 = omp_get_wtime();	
    for(int k = 0; k < REPS; k++)
    {   
        precondition(diag, sizes[0], initialResidual, reciprocalDiagPtr);
    }
    t1 = omp_get_wtime();
    cout << "Time precondition: " << (t1 - t0) / REPS << endl;
    
    cout << "Threads: " << omp_get_max_threads() << endl;
    
    // Check results
    reciprocalDiagPtrCopy = new double[sizes[0]];
    
    for (int cell=0; cell<sizes[0]; cell++)
    {
        reciprocalDiagPtrCopy[cell] = 1.0/diagCopy[cell];
    }
    
    for (int cell=0; cell<sizes[0]; cell++)
    {
        diagCopy[cell] = reciprocalDiagPtrCopy[cell] * initialResidual[cell];
    }
    
    for(int j = 0; j < sizes[0]; j++)
    {
        
        if(fabs(reciprocalDiagPtr[j] - reciprocalDiagPtrCopy[j]) > 1e-12)
        {
            cout << "Index: " << j << " error in reciprocalDiagPtr == " << reciprocalDiagPtr[j] << " and reciprocalDiagPtrCopy == " << reciprocalDiagPtrCopy[j] << endl;
            break;
        }
        if(fabs(diag[j] - diagCopy[j]) > 1e-12)
        {
            cout << "Index: " << j << " error in diag == " << diag[j] << " and diagCopy == " << diagCopy[j] << endl;
            break;
        }
    }

    //Destroy matrix
    destroy(sizes, psi, diag, upper, lower, upperAddr, lowerAddr, initialResidual, reciprocalDiagPtr, ownerStartAddr);

    return 0;
}