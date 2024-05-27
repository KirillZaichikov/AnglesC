#include <cmath>
#include <cstring>
#include <iostream>


double dotProduct(double* firstVec, double* secondVec, int length) {
    double res = 0;
    for (int i = 0; i < length; i++)
        res += firstVec[i] * secondVec[i];
    return res;
}

double vecNorm(double* vec, int length) {
    double norm = 0;
    for (int i = 0; i < length; i++)
        norm += vec[i] * vec[i];
    return sqrt(norm);
}

void normalizeVecs(double* vecs, int dimension, int numVec, double eps) {
    double norm = 0;
    for (int i = 0; i < numVec; i++) {
        norm = vecNorm(&vecs[i * dimension], dimension);
        for (int j = 0; j < dimension; j++)
            vecs[i * dimension + j] = (vecs[i * dimension + j] / norm) * eps;
    }
}

double getAngle(double* vec1, double* vec2, double dimension) {
    return acos((dotProduct(vec1, vec2, dimension)) / (vecNorm(vec1, dimension) * vecNorm(vec2, dimension)));
}

void ortVecs(double* vecs, int dimension, int numVecs, double* projSum) {
    double projCoeff = 0;

    for (int i = 1; i < numVecs; i++) {
        for (int j = 0; j < i; j++) {
            projCoeff = dotProduct(&vecs[j * dimension], &vecs[i * dimension], dimension) / dotProduct(&vecs[j * dimension], &vecs[j * dimension], dimension);
            for (int k = 0; k < dimension; k++)
                projSum[k] += projCoeff * vecs[j * dimension + k];
        }
        for (int j = 0; j < dimension; j++) {
            vecs[i * dimension + j] -= projSum[j];
            projSum[j] = 0;
        }
    }
}

void matrixTranspose(double* matr, int rows, int cols) {
    double temp;

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) {
            if (j > i) {
                temp = matr[i * cols + j];
                matr[i * cols + j] = matr[j * cols + i];
                matr[j * cols + i] = temp;
            }
        }
}

void matrixMult(double* res, double* A, int rowsA, int colsA, double* B, int rowsB, int colsB) {
    for (int i = 0; i < rowsA; i++)
        for (int j = 0; j < colsB; j++) {
            res[i * colsB + j] = 0;
            for (int k = 0; k < colsA; k++)
                res[i * colsB + j] += (A[i * colsA + k] * B[k * colsB + j]);
        }
}

double dotProductColMaj(double* firstVec, double* secondVec, int rows, int cols) {
    double res = 0;
    for (int i = 0; i < rows; i++)
        res += firstVec[i * cols] * secondVec[i * cols];
    return res;
}

double vecNormColMaj(double* vec, int rows, int cols) {
    double norm = 0;
    for (int i = 0; i < rows; i++)
        norm += vec[i * cols] * vec[i * cols];
    return sqrt(norm);
}

void normalizeVecsColMaj(double* vecs, int rows, int cols, double eps) {
    double norm = 0;
    for (int i = 0; i < rows; i++) {
        norm = vecNormColMaj(&vecs[i], rows, cols);
        for (int j = 0; j < cols; j++)
            vecs[i + j * cols] = (vecs[i + j * cols] / norm) * eps;
    }
}

void ortVecsColMaj(double* vecs, int rows, int cols, double* projSum) {
    double projCoeff = 0;

    for (int i = 1; i < cols; i++) {
        for (int j = 0; j < i; j++) {
            projCoeff = dotProductColMaj(&vecs[j], &vecs[i], rows, cols) / dotProductColMaj(&vecs[j], &vecs[j], rows, cols);
            for (int k = 0; k < rows; k++)
                projSum[k] += projCoeff * vecs[j + k * cols];
        }
        for (int j = 0; j < rows; j++) {
            vecs[i + j * cols] -= projSum[j];
            projSum[j] = 0;
        }
    }
}

void solver(double* M, int dim, double* Res){
    double D = pow((-M[3]-M[0]), 2.0) - 4 * (M[0]*M[3]-M[1]*M[2]);
    Res[0] = ((M[3]+M[0])+pow(D, 0.5))/2;
    Res[1] = ((M[3]+M[0])-pow(D, 0.5))/2;
    //std::cout << Res[0]<< ' ' << Res[1] << std::endl;
}

void qrDecomposition(double* R, double* Q, double* A, int rows, int cols, double* projSum) {
    memcpy(Q, A, cols * cols * sizeof(double));
    matrixTranspose(Q, cols, cols);
    ortVecs(Q, cols, rows, projSum);
    normalizeVecs(Q, cols, rows, 1.0);
    matrixMult(R, Q, cols, cols, A, cols, cols);
    matrixTranspose(Q, cols, cols);
}

void findEigenValues(double* A, double *R, double *Q, double* projSum, int rows, int cols, int numberOfIterations, double *eigenValues) {
    /*
    for (int32_t i = 0; i < rows; i++)
        for (int32_t j = 0; j < rows; j++) {
            printf("%lf ", A[i*cols + j]);
        }*/
    
    for (int i = 0; i < numberOfIterations; i++) {
        // printf("%d",i);
        qrDecomposition(R, Q, A, rows, cols, projSum);
        matrixMult(A, R, cols, cols, Q, cols, cols);
        /*
        for (int32_t i = 0; i < rows; i++)
            for (int32_t j = 0; j < rows; j++) {
                printf("%lf ", A[i*cols + j]);
            }*/
    }
    for (int i = 0; i < rows; i++) {
        eigenValues[i] = A[i * cols + i];
        //printf("%lf", eigenValues[i]);
    }
}