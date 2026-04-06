#include <iostream>
#include <cstdlib>
#include <cstdio>
using namespace std;

void init(
    int nCells,
    double *diagPtr,
    double **reciprocalDiagPtr,
    double *psi)
{
    *reciprocalDiagPtr = new double[nCells];

    for (int cell = 0; cell < nCells; cell++)
    {
        (*reciprocalDiagPtr)[cell] = 1.0 / diagPtr[cell];
    }
}

void precondition(
    double *diag,
    const int nCells,
    const int nFaces,
    const int *upperAddr,
    const int *lowerAddr,
    double *upper,
    double *lower,
    const double *initialResidual,
    const double *reciprocalDiagPtr)
{
    for (int cell = 0; cell < nCells; cell++)
    {
        diag[cell] = reciprocalDiagPtr[cell] * diag[cell];
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
    int **ownerStartAddr
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
    string filenameAfter = "../MatricesConXIT" + to_string(simpleFoamItDirectory) + "/matrixAfterDiag.bin";
    string filenameAfter2 = "../MatricesConXIT" + to_string(simpleFoamItDirectory) + "/matrixAfter2Diag.bin";
    string filenameAfter3 = "../MatricesConXIT" + to_string(simpleFoamItDirectory) + "/matrixAfter3Diag.bin";

    // Read matrix
    binaryRead(filenameBefore.c_str(), &sizes, &diag, &upperAddr, &lowerAddr, &upper, &lower, &initialResidual, &psi, &ownerStartAddr);

    // Init preconditioner
    init(sizes[0], diag, &reciprocalDiagPtr, psi);

    // Call preconditioner
    precondition(diag, sizes[0], sizes[1], upperAddr, lowerAddr, upper, lower, initialResidual, reciprocalDiagPtr);

    // Print matrix
    binaryPrint(filenameAfter.c_str(), sizes, diag, upperAddr, lowerAddr, upper, lower, initialResidual, psi, ownerStartAddr);

    // Call 2 time preconditioner
    precondition(diag, sizes[0], sizes[1], upperAddr, lowerAddr, upper, lower, initialResidual, reciprocalDiagPtr);

    // Print matrix
    binaryPrint(filenameAfter2.c_str(), sizes, diag, upperAddr, lowerAddr, upper, lower, initialResidual, psi, ownerStartAddr);

    // Call 3 time preconditioner
    precondition(diag, sizes[0], sizes[1], upperAddr, lowerAddr, upper, lower, initialResidual, reciprocalDiagPtr);

    // Print matrix
    binaryPrint(filenameAfter3.c_str(), sizes, diag, upperAddr, lowerAddr, upper, lower, initialResidual, psi, ownerStartAddr);

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