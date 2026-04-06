#include <iostream>
#include <cstdlib>
#include<cstdio>
using namespace std;

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
    for (int cell=0; cell<nCells; cell++)
    {
        (*reciprocalDiagPtr)[cell] = diagPtr[cell];
    }

    // Calculate the DIC diagonal
    for (int face=0; face<nFaces; face++)
    {
        (*reciprocalDiagPtr)[upperAddr[face]] -= upper[face] * upper[face] / (*reciprocalDiagPtr)[lowerAddr[face]];
    }


    // Calculate the reciprocal of the preconditioned diagonal
    for (int cell=0; cell<nCells; cell++)
    {
        (*reciprocalDiagPtr)[cell] = 1.0/(*reciprocalDiagPtr)[cell];
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

    for (int cell=0; cell<nCells; cell++)
    {
        diagPtr[cell] = reciprocalDiagPtr[cell]*initialResidual[cell];
    }

    for (int face=0; face<nFaces; face++)
    {
        diagPtr[upperAddr[face]] -= reciprocalDiagPtr[upperAddr[face]]*upper[face]*diagPtr[lowerAddr[face]];
    }

    for (int face=nFacesM1; face>=0; face--)
    {
        diagPtr[lowerAddr[face]] -= reciprocalDiagPtr[lowerAddr[face]]*upper[face]*diagPtr[upperAddr[face]];
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
    double** psi
)
{
            FILE *binaryFile = fopen(fileName, "rb");   
            
            *sizes = new int[7];
            fread(*sizes, sizeof(int), 7, binaryFile);
            *diag = new double[(*sizes)[0]];
            *lower = new double[(*sizes)[1]];
            *lowerAddr = new int[(*sizes)[2]];
            *upper = new double[(*sizes)[3]];
            *upperAddr = new int[(*sizes)[4]];
            *initialResidual = new double[(*sizes)[5]];
            *psi = new double[(*sizes)[6]];

            fread(*diag           , sizeof(double), (*sizes)[0] , binaryFile);
            fread(*lower          , sizeof(double), (*sizes)[1] , binaryFile);
            fread(*lowerAddr      , sizeof(int)   , (*sizes)[2] , binaryFile);
            fread(*upper          , sizeof(double), (*sizes)[3] , binaryFile);
            fread(*upperAddr      , sizeof(int)   , (*sizes)[4] , binaryFile);
            fread(*initialResidual, sizeof(double)   , (*sizes)[5] , binaryFile);
            fread(*psi            , sizeof(double)   , (*sizes)[6] , binaryFile);
          
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
    const double* psi
)
{
            FILE *binaryFile = fopen(fileName, "wb");   

            fwrite(sizes          , sizeof(int)   , 7       , binaryFile);
            fwrite(diag           , sizeof(double), sizes[0], binaryFile);
            fwrite(lower          , sizeof(double), sizes[1], binaryFile);
            fwrite(lowerAddr      , sizeof(int)   , sizes[2], binaryFile);
            fwrite(upper          , sizeof(double), sizes[3], binaryFile);
            fwrite(upperAddr      , sizeof(int)   , sizes[4], binaryFile);
            fwrite(initialResidual, sizeof(double), sizes[5], binaryFile);
            fwrite(psi            , sizeof(double), sizes[6], binaryFile);
          
            fclose(binaryFile);
}

void textPrint
(
    const char* fileName, 
    const int* sizes,
    const double* diag,
    const int* upperAddr,
    const int* lowerAddr,
    const double* upper,
    const double* lower,
    const double* initialResidual,
    const double* psi
)
{
            FILE *textFile = fopen(fileName, "w");

            for(int i = 0; i < sizes[0]; i++)
            {
                fprintf(textFile, "%f ", diag[i]);
            }
            fprintf(textFile, "\n");
          
            for(int i = 0; i < sizes[1]; i++)
            {
                fprintf(textFile, "%f ", lower[i]);
            }
            fprintf(textFile, "\n");
          
            for(int i = 0; i < sizes[2]; i++)
            {
                fprintf(textFile, "%d ", lowerAddr[i]);
            }
            fprintf(textFile, "\n");
          
            for(int i = 0; i < sizes[3]; i++)
            {
                fprintf(textFile, "%f ", upper[i]);
            }
            fprintf(textFile, "\n");
          
            for(int i = 0; i < sizes[4]; i++)
            {
                fprintf(textFile, "%d ", upperAddr[i]);
            }
            fprintf(textFile, "\n");
          
            for(int i = 0; i < sizes[5]; i++)
            {
                fprintf(textFile, "%f ", initialResidual[i]);
            }
            fprintf(textFile, "\n");
          
            for(int i = 0; i < sizes[6]; i++)
            {
                fprintf(textFile, "%f ", psi[i]);
            }
            fprintf(textFile, "\n");
          
            fclose(textFile);
}

int main(){
    int *sizes;
    double *diag;
    double *lower;
    double *upper;
    int *lowerAddr;
    int *upperAddr;
    double *reciprocalDiagPtr;
    double *initialResidual;
    double *psi;

    //Read matrix values
    binaryRead("../MatricesConXIT500/matrixBefore.bin", &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi);

    //Init algorithm
    init(sizes[0], sizes[3], upper, upperAddr, lowerAddr, diag, &reciprocalDiagPtr);

    // Precondition calling
    precondition(diag, sizes[0], sizes[3], upper, upperAddr, lowerAddr, initialResidual, reciprocalDiagPtr);

    //Print matrix
    //textPrint("../MatricesConX/matrixAfterDiag.txt", sizes, diag, upperAddr, lowerAddr, upper, lower, initialResidual, psi); 
    binaryPrint("../MatricesConXIT500/matrixAfterDIC.bin", sizes, diag, upperAddr, lowerAddr, upper, lower, initialResidual, psi); 

    //Destroy matrix
    destroy(sizes, psi, diag, upper, lower, upperAddr, lowerAddr, initialResidual, reciprocalDiagPtr);

    return 0;
}
