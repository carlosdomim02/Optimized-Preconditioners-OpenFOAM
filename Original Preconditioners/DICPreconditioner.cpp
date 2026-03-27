#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
using namespace std;

void init(
    int nCells,
    int nFaces,
    double *upper,
    int *upperAddr,
    int *lowerAddr,
    double *diagPtr,
    double **reciprocalDiagPtr)
{
    (*reciprocalDiagPtr) = new double[nCells];

    for (int cell = 0; cell < nCells; cell++)
    {
        (*reciprocalDiagPtr)[cell] = diagPtr[cell];
    }

    for (int face = 0; face < nFaces; face++)
    {
        (*reciprocalDiagPtr)[upperAddr[face]] -= upper[face] * upper[face] / (*reciprocalDiagPtr)[lowerAddr[face]];
    }

    for (int cell = 0; cell < nCells; cell++)
    {
        (*reciprocalDiagPtr)[cell] = 1.0 / (*reciprocalDiagPtr)[cell];
    }
}

void precondition(
    double *diagPtr,
    const int nCells,
    const int nFaces,
    const double *upper,
    const int *upperAddr,
    const int *lowerAddr,
    const double *initialResidual,
    const double *reciprocalDiagPtr)
{
    int nFacesM1 = nFaces - 1;

    for (int cell = 0; cell < nCells; cell++)
    {
        diagPtr[cell] = reciprocalDiagPtr[cell] * diagPtr[cell];
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
    const int *sizes,
    const double *psi,
    const double *diag,
    const double *upper,
    const double *lower,
    const int *upperAddr,
    const int *lowerAddr,
    const double *initialResidual,
    const double *reciprocalDiagPtr,
    const int *ownerStartAddr)
{
    delete[] (psi);
    delete[] (diag);
    delete[] (sizes);
    delete[] (upper);
    delete[] (lower);
    delete[] (upperAddr);
    delete[] (lowerAddr);
    delete[] (ownerStartAddr);
    delete[] (initialResidual);
    delete[] (reciprocalDiagPtr);
}

void binaryRead(
    const char *fileName,
    int **sizes,
    double **diag,
    int **upperAddr,
    int **lowerAddr,
    double **upper,
    double **lower,
    double **initialResidual,
    double **psi,
    int **ownerStartAddr)
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

void binaryPrint(
    const char *fileName,
    const int *sizes,
    const double *diag,
    const int *upperAddr,
    const int *lowerAddr,
    const double *upper,
    const double *lower,
    const double *initialResidual,
    const double *psi,
    const int *ownerStartAddr)
{
    FILE *binaryFile = fopen(fileName, "wb");

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

void execInitAndPreconditionAt(const int simpleFoamItDirectory)
{
    int *sizes;
    double *diag;
    double *lower;
    double *upper;
    int *lowerAddr;
    int *upperAddr;
    int *ownerStartAddr;
    double *reciprocalDiagPtr;
    double *initialResidual;
    double *psi;

    string filenameBefore = "../MatricesConXIT" + to_string(simpleFoamItDirectory) + "/matrixBefore.bin";
    string filenameAfter = "../MatricesConXIT" + to_string(simpleFoamItDirectory) + "/matrixAfterDIC.bin";
    // Read matrix values
    binaryRead(filenameBefore.c_str(), &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi, &ownerStartAddr);

    // Init
    init(sizes[0], sizes[3], upper, upperAddr, lowerAddr, diag, &reciprocalDiagPtr);

    // Precondition
    precondition(diag, sizes[0], sizes[3], upper, upperAddr, lowerAddr, initialResidual, reciprocalDiagPtr);

    // Print matrix
    binaryPrint(filenameAfter.c_str(), sizes, diag, upperAddr, lowerAddr, upper, lower, initialResidual, psi, ownerStartAddr);

    // Destroy matrix
    destroy(sizes, psi, diag, upper, lower, upperAddr, lowerAddr, initialResidual, reciprocalDiagPtr, ownerStartAddr);
}

int main()
{
    execInitAndPreconditionAt(250);
    for (int i = 100; i <= 500; i += 100)
    {
        execInitAndPreconditionAt(i);
    } 

    return 0;
}
