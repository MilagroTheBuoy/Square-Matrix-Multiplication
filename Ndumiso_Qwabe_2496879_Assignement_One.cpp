#include <bits/stdc++.h> 
#include <iostream>
#include <vector>

void printMatrix(std::vector<std::vector<float>> matrix, size_t n){
    for(size_t i = 0; i < n; i++){

        for(size_t j = 0; j < n; j++){

            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

//eg. 10 is a power of two, 10 & 01 = 00. negated is 11
bool isPowerOfTwo(size_t num){
    return !(num & (num - 1));
}

void inflateMatrix(std::vector<std::vector<float>> &matrix){
    size_t oldRowSize = matrix.size();
    size_t oldColSize = matrix[0].size();
    size_t newRowExponent, newColExponent;

    newRowExponent = std::ceil(log2(oldRowSize));
    newColExponent = std::ceil(log2(oldColSize));

    if(newRowExponent > newColExponent){
        newColExponent = newRowExponent;
    }
    else if(newRowExponent < newColExponent){
        newRowExponent = newColExponent;
    }

    size_t newRowSize = pow(2, newRowExponent);
    size_t newColSize = pow(2, newColExponent);
    size_t rowsToAdd = newRowSize - oldRowSize;
    size_t colsToAdd = newColSize - oldColSize;

    for (int i = 0; i < oldRowSize; i++){
        for (int j = 0; j < colsToAdd; j++){
            matrix[i].push_back(0);
        }
    }

    for (int k = 0; k < rowsToAdd; k++){
        matrix.push_back(std::vector<float>(newRowSize));
    }
}

void deflateMatrix(std::vector<std::vector<float>> &matrix, size_t n){
    size_t currSize = matrix.size();
    size_t amountToDelete = currSize - n;

    for(size_t i = 0; i < amountToDelete; i++){
        matrix.pop_back();
    }

    for(size_t i = 0; i < n; i++){
        for(size_t j = 0; j < amountToDelete; j++){
            matrix[i].pop_back();
        }
    }
}

std::vector<std::vector<float>> square_matrix_mult(std::vector<std::vector<float>> 
                                                &A, std::vector<std::vector<float>> &B){

    size_t n = A.size();
    std::vector<std::vector<float>> C(n, std::vector<float>(n));

    for(size_t i = 0; i < n; i++){

        for(size_t j = 0; j < n; j++){

            for(size_t k = 0; k < n; k++){
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return C;
}

std::vector<std::vector<float>> addMatrices(const std::vector<std::vector<float>> 
                                        &A, const std::vector<std::vector<float>> &B){

    size_t matrix_size = A.size();
    std::vector<std::vector<float>> C(matrix_size, std::vector<float>(matrix_size));

    for (size_t i = 0; i < matrix_size; i++){
        for (size_t j = 0; j < matrix_size; j++){
            C[i][j] += A[i][j] + B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<float>> subtractMatrices(const std::vector<std::vector<float>> 
                                        &A, const std::vector<std::vector<float>> &B){

    size_t matrix_size = A.size();
    std::vector<std::vector<float>> C(matrix_size, std::vector<float>(matrix_size));

    for (size_t i = 0; i < matrix_size; i++){
        for (size_t j = 0; j < matrix_size; j++){
            C[i][j] += A[i][j] - B[i][j];
        }
    }
    return C;
}

// 0 is rowStart, 1 is rowEnd, 2 is colStart, 3 is colEnd
/*std::vector<std::vector<float>> addMatrices(const std::vector<std::vector<float>> 
                                        &A, const std::vector<std::vector<float>> &B,
                                        int A_indices[], int B_indices[], int matrix_size){

    std::vector<std::vector<float>> C(matrix_size, std::vector<float>(matrix_size));

    if(A_indices[0] != -1){
        for (size_t i = A_indices[0]; i < A_indices[1] + 1; i++){
            for (size_t j = A_indices[2]; j < A_indices[3] + 1; j++){
                C[i - A_indices[0]][j - A_indices[2]] += A[i][j];
            }
        }
    }
    else{
        for (size_t i = 0; i < C.size(); i++){
            for (size_t j = 0; j < C.size(); j++){
                C[i][j] += A[i][j];
            }
        }
    }

    if(B_indices[0] != -1){
        for (size_t i = B_indices[0]; i < B_indices[1] + 1; i++){
            for (size_t j = B_indices[2]; j < B_indices[3] + 1; j++){
                C[i - B_indices[0]][j - B_indices[2]] += B[i][j];
            }
        }
    }
    else{
        for (size_t i = 0; i < C.size(); i++){
            for (size_t j = 0; j < C.size(); j++){
                C[i][j] += B[i][j];
            }
        }
    }
    return C;
}

std::vector<std::vector<float>> subtractMatrices(const std::vector<std::vector<float>> 
                                        &A, const std::vector<std::vector<float>> &B,
                                        int A_indices[], int B_indices[], int matrix_size){

    std::vector<std::vector<float>> C(matrix_size, std::vector<float>(matrix_size));

    if(A_indices[0] != -1){
        for (size_t i = A_indices[0]; i < A_indices[1] + 1; i++){
            for (size_t j = A_indices[2]; j < A_indices[3] + 1; j++){
                //int d = i - A_indices[0];
                //int g = j - A_indices[2];
                C[i - A_indices[0]][j - A_indices[2]] += A[i][j];
            }
        }
    }
    else{
        for (size_t i = 0; i < C.size(); i++){
            for (size_t j = 0; j < C.size(); j++){
                C[i][j] += A[i][j];
            }
        }
    }

    if(B_indices[0] != -1){
        for (size_t i = B_indices[0]; i < B_indices[1] + 1; i++){
            for (size_t j = B_indices[2]; j < B_indices[3] + 1; j++){
                C[i - B_indices[0]][j - B_indices[2]] -= B[i][j];
            }
        }
    }
    else{
        for (size_t i = 0; i < C.size(); i++){
            for (size_t j = 0; j < C.size(); j++){
                C[i][j] -= B[i][j];
            }
        }
    }

    return C;
}*/

std::vector<std::vector<float>> combineIntoC(std::vector<std::vector<float>> &C, 
                std::vector<std::vector<float>> C_00, std::vector<std::vector<float>> C_01, 
                std::vector<std::vector<float>> C_10, std::vector<std::vector<float>> C_11){

    size_t subMatrixSize = C_00.size();
    for(int i = 0; i < subMatrixSize; i++){
        C_00[i].insert(C_00[i].end(), C_01[i].begin(), C_01[i].end());
        C_10[i].insert(C_10[i].end(), C_11[i].begin(), C_11[i].end());
    }
    C_00.insert(C_00.end(), C_10.begin(), C_10.end());
    return C_00;
}

// 0 is rowStart, 1 is rowEnd, 2 is colStart, 3 is colEnd
std::vector<std::vector<float>> square_matrix_mult_recurse(const std::vector<std::vector<float>> 
                                                            &A, const std::vector<std::vector<float>> &B, 
                                                            int A_indices[], int B_indices[]){
    
    size_t n = (A_indices[1] - A_indices[0]) + 1;
    std::vector<std::vector<float>> C(n, std::vector<float>(n));
    if(n == 1){
        size_t rowA = A_indices[0];
        size_t colA = A_indices[2];
        size_t rowB = B_indices[0];
        size_t colB = B_indices[2];
        C[0][0] = A[rowA][colA]*B[rowB][colB];
    }
    else{
        n = n/2;

        int A_00[4] = {A_indices[0], A_indices[0] + n - 1, A_indices[2], A_indices[2] + n - 1};
        int A_01[4] = {A_indices[0], A_indices[0] + n - 1, A_indices[2] + n, A_indices[3]};
        int A_10[4] = {A_indices[0] + n, A_indices[1], A_indices[2], A_indices[2] + n -1};
        int A_11[4] = {A_indices[0] + n, A_indices[1], A_indices[2] + n, A_indices[3]};
        int B_00[4] = {B_indices[0], B_indices[0] + n - 1, B_indices[2], B_indices[2] + n - 1};
        int B_01[4] = {B_indices[0], B_indices[0] + n - 1, B_indices[2] + n, B_indices[3]};
        int B_10[4] = {B_indices[0] + n, B_indices[1], B_indices[2], B_indices[2] + n -1};
        int B_11[4] = {B_indices[0] + n, B_indices[1], B_indices[2] + n, B_indices[3]};

        std::vector<std::vector<float>> C_00 = addMatrices(square_matrix_mult_recurse(A, B, A_00, B_00), 
                                                square_matrix_mult_recurse(A, B, A_01, B_10));
        std::vector<std::vector<float>> C_01 = addMatrices(square_matrix_mult_recurse(A, B, A_00, B_01), 
                                                square_matrix_mult_recurse(A, B, A_01, B_11));
        std::vector<std::vector<float>> C_10 = addMatrices(square_matrix_mult_recurse(A, B, A_10, B_00), 
                                                square_matrix_mult_recurse(A, B, A_11, B_10));
        std::vector<std::vector<float>> C_11 = addMatrices(square_matrix_mult_recurse(A, B, A_10, B_01), 
                                                square_matrix_mult_recurse(A, B, A_11, B_11));

        C = combineIntoC(C, C_00, C_01, C_10, C_11);
    }
    return C;
}

std::vector<std::vector<float>> createSubMatrix(const std::vector<std::vector<float>> &A, int indices[], int subArraySize){
    
    std::vector<std::vector<float>> subArray(subArraySize, std::vector<float>(subArraySize));

    for (size_t i = indices[0]; i < indices[1] + 1; i++){
        for (size_t j = indices[2]; j < indices[3] + 1; j++){
            subArray[i - indices[0]][j - indices[2]] += A[i][j];
        }
    }
    return subArray;
}

// 0 is rowStart, 1 is rowEnd, 2 is colStart, 3 is colEnd
std::vector<std::vector<float>> square_matrix_mult_strassen(std::vector<std::vector<float>> 
                                                            &A, std::vector<std::vector<float>> &B){
    size_t n = A.size();
    std::vector<std::vector<float>> C(n, std::vector<float>(n));
    if(n == 1){
        //size_t rowA = A_indices[0];
        //size_t colA = A_indices[2];
        //size_t rowB = B_indices[0];
        //size_t colB = B_indices[2];
        C[0][0] = A[0][0]*B[0][0];
    }
    else{
        n /= 2;
        int ignore[4] = {-1, -1, -1, -1};
        int A_00[4] = {0, n - 1, 0, n - 1};
        int A_01[4] = {0, n - 1, n, 2*n - 1};
        int A_10[4] = {n, 2*n -1, 0, n -1};
        int A_11[4] = {n, 2*n - 1, n, 2*n - 1};
        int B_00[4] = {0, n - 1, 0, n - 1};
        int B_01[4] = {0, n - 1, n, 2*n - 1};
        int B_10[4] = {n, 2*n -1, 0, n -1};
        int B_11[4] = {n, 2*n - 1, n, 2*n - 1};
        
        std::vector<std::vector<float>> A_00_mat = createSubMatrix(A, A_00, n);
        std::vector<std::vector<float>> A_01_mat = createSubMatrix(A, A_01, n);
        std::vector<std::vector<float>> A_10_mat = createSubMatrix(A, A_10, n);
        std::vector<std::vector<float>> A_11_mat = createSubMatrix(A, A_11, n);
        std::vector<std::vector<float>> B_00_mat = createSubMatrix(B, B_00, n);
        std::vector<std::vector<float>> B_01_mat = createSubMatrix(B, B_01, n);
        std::vector<std::vector<float>> B_10_mat = createSubMatrix(B, B_10, n);
        std::vector<std::vector<float>> B_11_mat = createSubMatrix(B, B_11, n);

        std::vector<std::vector<float>> S1 = subtractMatrices(B_01_mat, B_11_mat);
        std::vector<std::vector<float>> S2 = addMatrices(A_00_mat, A_01_mat);
        std::vector<std::vector<float>> S3 = addMatrices(A_10_mat, A_11_mat);
        std::vector<std::vector<float>> S4 = subtractMatrices(B_10_mat, B_00_mat);
        std::vector<std::vector<float>> S5 = addMatrices(A_00_mat, A_11_mat);
        std::vector<std::vector<float>> S6 = addMatrices(B_00_mat, B_11_mat);
        std::vector<std::vector<float>> S7 = subtractMatrices(A_01_mat, A_11_mat);
        std::vector<std::vector<float>> S8 = addMatrices(B_10_mat, B_11_mat);
        std::vector<std::vector<float>> S9 = subtractMatrices(A_00_mat, A_10_mat);
        std::vector<std::vector<float>> S10 = addMatrices(B_00_mat, B_01_mat);

        std::vector<std::vector<float>> P1 = square_matrix_mult_strassen(A_00_mat, S1);
        std::vector<std::vector<float>> P2 = square_matrix_mult_strassen(S2, B_11_mat);
        std::vector<std::vector<float>> P3 = square_matrix_mult_strassen(S3, B_00_mat);
        std::vector<std::vector<float>> P4 = square_matrix_mult_strassen(A_11_mat, S4);
        std::vector<std::vector<float>> P5 = square_matrix_mult_strassen(S5, S6);
        std::vector<std::vector<float>> P6 = square_matrix_mult_strassen(S7, S8);
        std::vector<std::vector<float>> P7 = square_matrix_mult_strassen(S9, S10);

        std::vector<std::vector<float>> C_00 = addMatrices(P5, P4);
        C_00 = subtractMatrices(C_00, P2);
        C_00 = addMatrices(C_00, P6);

        std::vector<std::vector<float>> C_01 = addMatrices(P1, P2);

        std::vector<std::vector<float>> C_10 = addMatrices(P3, P4);

        std::vector<std::vector<float>> C_11 = addMatrices(P5, P1);
        C_11 = subtractMatrices(C_11, P3);
        C_11 = subtractMatrices(C_11, P7);

        C = combineIntoC(C, C_00, C_01, C_10, C_11);

    }
    return C;
}


int main(){

    size_t originalSize;
    std::vector<std::vector<float>> A = {{1, 5, 9}, {6, 4, 8}, {6, 4, 7}};
    std::vector<std::vector<float>> B = {{4, 5, 9}, {8, 2, 6}, {6, 4, 5}};

    std::vector<std::vector<float>> C = square_matrix_mult(A,B);
    
    originalSize = A.size();
    inflateMatrix(A);
    inflateMatrix(B);

    int initialA_Indices[] = {0, A.size() - 1, 0, A.size()- 1};
    int initialB_Indices[] = {0, B.size() - 1, 0, B.size() - 1};

    std::vector<std::vector<float>> D = square_matrix_mult_recurse(A,B, initialA_Indices, initialB_Indices);
    std::vector<std::vector<float>> E = square_matrix_mult_strassen(A,B);
    
    deflateMatrix(C, originalSize);
    deflateMatrix(D, originalSize);
    deflateMatrix(E, originalSize);

    printMatrix(C, C.size());
    printMatrix(D, D.size());
    printMatrix(E, E.size());
    /* 44 15 
        56 38 */
    return 0;
}