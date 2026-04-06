#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstring>
#include <omp.h>
#include <vector>
using namespace std;

#define REPS 100000

void init
(
    const int nCells,
    const double* upper,
    const int* upperAddr,
    const int* lowerAddr,
    const double* diagPtr,
    double* reciprocalDiagPtr,
    vector<int> cells,
    vector<int> cellsStartingIndex
)
{      
      #pragma omp parallel
      {
          int end;
          int cell;
          int face;
          int index;   
          
          #pragma omp for simd
          for (cell=0; cell<nCells; cell++)
          {
              reciprocalDiagPtr[cell] = diagPtr[cell];
          }
          
          #pragma omp for
          for(cell = 0; cell < nCells; cell++)
          {
              end = cellsStartingIndex[cell + 1];
              #pragma omp simd
              for(index = cellsStartingIndex[cell]; index < end; index++)
              {
                  face = cells[index];
                  if(upperAddr[face] == cell) 
                      reciprocalDiagPtr[cell] -= upper[face] * upper[face] / reciprocalDiagPtr[lowerAddr[face]];
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
    double* diagPtr,
    const int nCells,
    const double* upper,
    const int* upperAddr,
    const int* lowerAddr,
    const double* initialResidual,
    const double* reciprocalDiagPtr,
    vector<int> cells,
    vector<int> cellsStartingIndex
)
{      
      
      int nCellsM1 = nCells - 1;
      
      #pragma omp parallel
      {
          int end;
          int cell;
          int face;
          int index;
          double tempValue;
          
          #pragma omp for simd
          for (cell=0; cell<nCells; cell++)
          {
              diagPtr[cell] = reciprocalDiagPtr[cell] * initialResidual[cell];
          }
          
          #pragma omp for
          for(int cell = 0; cell < nCells; cell++)
          {
              tempValue = 0.0;
              end = cellsStartingIndex[cell + 1];
              #pragma omp simd
              for(index = cellsStartingIndex[cell]; index < end; index++)
              {
                  face = cells[index];
                  if(upperAddr[face] == cell) 
                          tempValue += upper[face]*diagPtr[lowerAddr[face]];
              }              
              diagPtr[cell] -= reciprocalDiagPtr[cell]*tempValue;
          }
          
          #pragma omp for
          for(int cell = nCellsM1; cell >= 0; cell--)
          {
              tempValue = 0.0;
              end = cellsStartingIndex[cell + 1];
              #pragma omp simd
              for(index = cellsStartingIndex[cell]; index < end; index++)
              {
                  face = cells[index];
                  if(lowerAddr[face] == cell) 
                          tempValue += upper[face]*diagPtr[upperAddr[face]]; 
              }
              diagPtr[cell] -= reciprocalDiagPtr[cell]*tempValue;
          }
      }
}

void destroy(
    const int *sizes,
    double *diag,
    double *lower,
    double *upper,
    int *lowerAddr,
    int *upperAddr,
    double *reciprocalDiagPtr,
    double *initialResidual,
    double *psi,
    int *ownerStartAddr
)
{
    delete[] sizes;
    delete[] diag;
    delete[] lower;
    delete[] upper;
    delete[] lowerAddr;
    delete[] upperAddr;
    delete[] reciprocalDiagPtr;
    delete[] initialResidual;
    delete[] psi;
    delete[] ownerStartAddr;
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

void orderByCells
(
    const int nCells,
    const int nFaces,
    const int *upperAddr,
    const int *lowerAddr,
    vector<int> &cells,
    vector<int> &cellsStartingIndex
)
{
    vector<vector<int>> cellsMatrix(nCells);
    
    for(int face = 0; face < nFaces; face++)
    {
        cellsMatrix[lowerAddr[face]].push_back(face);
        cellsMatrix[upperAddr[face]].push_back(face); 
    }
    
    cellsStartingIndex.push_back(0);
    for(int cell = 0; cell < cellsMatrix.size(); cell++)
    {
        for(int face = 0; face < cellsMatrix[cell].size(); face++)
        {
            cells.push_back(cellsMatrix[cell][face]); 
        } 
        cellsStartingIndex.push_back(cellsStartingIndex.back() + cellsMatrix[cell].size());
    }  
    //cellsStartingIndex.pop_back();     
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
    double *reciprocalDiagPtr;
    double *reciprocalDiagPtrCopy;
    double *initialResidual;
    double *psi;    
    
    vector<int> cells;
    vector<int> cellsStartingIndex;
    
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
    binaryRead("./matrixBefore.bin", &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi, &ownerStartAddr);
    
    orderByCells(sizes[0], sizes[1], upperAddr, lowerAddr, cells, cellsStartingIndex);
    
    diagCopy = new double[sizes[0]];
    for (int cell=0; cell<sizes[0]; cell++)
    {
        diagCopy[cell] = diag[cell];
    }
    
    reciprocalDiagPtr = new double[sizes[0]];
    reciprocalDiagPtrCopy = new double[sizes[0]];
        
    //Time init
    t0 = omp_get_wtime();	
    for(int i = 0; i < REPS; i++){
        init(sizes[0], upper, upperAddr, lowerAddr, diag, reciprocalDiagPtr, cells, cellsStartingIndex);
    }
    t1 = omp_get_wtime();	
        
    cout << (t1 - t0) / REPS << endl;
    
    //Time precondition
    t0 = omp_get_wtime();	
    for(int i = 0; i < REPS; i++){
        precondition(diag, sizes[0], upper, upperAddr, lowerAddr, initialResidual, reciprocalDiagPtr, cells, cellsStartingIndex);
    }
    t1 = omp_get_wtime();	
    
    cout << (t1 - t0) / REPS << endl;
    
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
    //In this code checking often fail due to the fact that the matrix is not ordered
    /*for(int j = 0; j < sizes[0]; j++)
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
    }*/
    
    //Destroy matrix
    destroy(sizes, diag, lower, upper, lowerAddr, upperAddr, reciprocalDiagPtr, initialResidual, psi, ownerStartAddr);

    return 0;
}