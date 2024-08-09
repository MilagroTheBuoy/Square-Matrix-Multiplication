#include <bits/stdc++.h> 
#include <iostream>
#include <chrono>
#include <fstream>
#include <string>

const int TRIALS_PER_N_VALUE = 3, MIN_RECURSION_DEPTH = 0;

class SpecialArray{
    public:
        float ** rows;
        float *arrayData;
        int arrSize;
        bool defaultInitialisation = true;
        std::string name = "";
        
        SpecialArray(int n){
            arrSize = n;
            rows = new float*[n];
            arrayData = new float[n*n]{0};
            for(int i = 0; i < n; i++){
                rows[i] = arrayData + i * n;
            }
        }

        SpecialArray(int n, float *matrix){
            arrSize = n;
            rows = new float*[n];
            arrayData = matrix;
            for(int i = 0; i < n; i++){
                rows[i] = (arrayData + i * n);
            }
            defaultInitialisation = false;
        }

        ~SpecialArray(){
            if(!defaultInitialisation){
                bool f =(*rows == nullptr);
                delete[] arrayData;
                delete[] rows;
                arrayData = nullptr;
                rows = nullptr;
            }
            else{
                delete[] arrayData;
                delete[] rows;
                arrayData = nullptr;
                rows = nullptr;
            }
        }
};

void printMatrix(float *matrix, size_t n){
    for(size_t row = 0; row < n; row++){
        for(size_t col = 0; col < n; col++){
            std::cout << *((matrix + row * n) + col) << " ";
        }
        std::cout << std::endl << std::endl;
    }
}

void printMatrix(SpecialArray *matrix, size_t n){
    for(size_t row = 0; row < n; row++){
        for(size_t col = 0; col < n; col++){
            std::cout << matrix->rows[row][col] << " ";
        }
        std::cout << std::endl << std::endl;
    }
}

//eg. 10 is a power of two, 10 & 01 = 00. negated is 11
bool isPowerOfTwo(size_t num){
    return !(num & (num - 1));
}

void addMatrices(SpecialArray *A, SpecialArray *B, SpecialArray &C, 
                int mat_C[], int n){

    int startRow = mat_C[0], endRow = mat_C[1], startCol = mat_C[2], endCol = mat_C[3]; 
    for (size_t row = 0; row < n; row++){
        for (size_t col = 0; col < n; col++){
            C.rows[startRow + row][startCol + col] += A->rows[row][col] + B->rows[row][col];
        }
    }
}

void addMatrices(SpecialArray *A, SpecialArray *B, SpecialArray &C, int mat_A[], int mat_B[],
                int mat_C[], int n){

    int startRowC = mat_C[0], endRowC = mat_C[1], startColC = mat_C[2], endColC = mat_C[3];
    int startRowB = mat_B[0], endRowB = mat_B[1], startColB = mat_B[2], endColB = mat_B[3];
    int startRowA = mat_A[0], endRowA = mat_A[1], startColA = mat_A[2], endColA = mat_A[3];

    for (size_t row = 0; row < n; row++){
        for (size_t col = 0; col < n; col++){
            C.rows[startRowC + row][startColC + col] += A->rows[startRowA + row][startColA + col]
                                                        + B->rows[startRowB + row][startColB + col];
        }
    }
}


void subtractMatrices(SpecialArray *A, SpecialArray *B, SpecialArray &C, int mat_A[], int mat_B[],
                int mat_C[], int n){

    int startRowC = mat_C[0], endRowC = mat_C[1], startColC = mat_C[2], endColC = mat_C[3];
    int startRowB = mat_B[0], endRowB = mat_B[1], startColB = mat_B[2], endColB = mat_B[3];
    int startRowA = mat_A[0], endRowA = mat_A[1], startColA = mat_A[2], endColA = mat_A[3];

    for (size_t row = 0; row < n; row++){
        for (size_t col = 0; col < n; col++){
            C.rows[startRowC + row][startColC + col] += A->rows[startRowA + row][startColA + col]
                                                        - B->rows[startRowB + row][startColB + col];
        }
    }
}

void addSingleMatrix(SpecialArray *A, SpecialArray &C, int mat_A[], int mat_C[], int n){

    int startRowC = mat_C[0], endRowC = mat_C[1], startColC = mat_C[2], endColC = mat_C[3];
    int startRowA = mat_A[0], endRowA = mat_A[1], startColA = mat_A[2], endColA = mat_A[3];

    for (size_t row = 0; row < n; row++){
        for (size_t col = 0; col < n; col++){
            C.rows[startRowC + row][startColC + col] += A->rows[startRowA + row][startColA + col];
        }
    }
}


void subtractSingleMatrix(SpecialArray *A, SpecialArray &C, int mat_A[], int mat_C[], int n){

    int startRowC = mat_C[0], endRowC = mat_C[1], startColC = mat_C[2], endColC = mat_C[3];
    int startRowA = mat_A[0], endRowA = mat_A[1], startColA = mat_A[2], endColA = mat_A[3];

    for (size_t row = 0; row < n; row++){
        for (size_t col = 0; col < n; col++){
            C.rows[startRowC + row][startColC + col] -= A->rows[startRowA + row][startColA + col];
        }
    }
}

void square_matrix_mult(float *A, float *B, float *C, int n){
    for(size_t row = 0; row < n; row++){
        for(size_t col = 0; col < n; col++){
            for(size_t k = 0; k < n; k++){
                int a_num = *((A + row * n) + k);
                int b_num = *((B + k * n) + col);
                *((C + row * n) + col) +=  a_num * b_num;
            }
        }
    }
}

SpecialArray square_matrix_mult(SpecialArray *A, SpecialArray *B, int n){

    SpecialArray C = SpecialArray(n);
    for(size_t row = 0; row < n; row++){
        for(size_t col = 0; col < n; col++){
            for(size_t k = 0; k < n; k++){
                //int a_num = *((A + row * n) + k);
                //int b_num = *((B + k * n) + col);
                C.rows[row][col] +=  A->rows[row][k] * B->rows[k][col];
            }
        }
    }
    return C;
}

// 0 is rowStart, 1 is rowEnd, 2 is colStart, 3 is colEnd
SpecialArray square_matrix_mult_recurse(SpecialArray *A, SpecialArray *B, int A_indices[], 
                                        int B_indices[], int arraySize){
    
    SpecialArray C = SpecialArray(arraySize);

    if(arraySize == 1){
        size_t rowA = A_indices[0];
        size_t colA = A_indices[2];
        size_t rowB = B_indices[0];
        size_t colB = B_indices[2];
        C.rows[0][0] = A->rows[rowA][colA] * B->rows[rowB][colB];
    }
    else{
        arraySize /= 2;

        int A_00[4] = {A_indices[0], A_indices[0] + arraySize - 1, A_indices[2],
                        A_indices[2] + arraySize - 1};
        int A_01[4] = {A_indices[0], A_indices[0] + arraySize - 1, A_indices[2] + arraySize, 
                        A_indices[3]};
        int A_10[4] = {A_indices[0] + arraySize, A_indices[1], A_indices[2], 
                        A_indices[2] + arraySize - 1};
        int A_11[4] = {A_indices[0] + arraySize, A_indices[1], A_indices[2] + arraySize,
                        A_indices[3]};
        int B_00[4] = {B_indices[0], B_indices[0] + arraySize - 1, B_indices[2],
                        B_indices[2] + arraySize - 1};
        int B_01[4] = {B_indices[0], B_indices[0] + arraySize - 1, B_indices[2] + arraySize, 
                        B_indices[3]};
        int B_10[4] = {B_indices[0] + arraySize, B_indices[1], B_indices[2], 
                        B_indices[2] + arraySize - 1};
        int B_11[4] = {B_indices[0] + arraySize, B_indices[1], B_indices[2] + arraySize,
                        B_indices[3]};

        int C_00[4] = {0, arraySize - 1, 0, arraySize - 1};
        int C_01[4] = {0, arraySize - 1, arraySize, 2*arraySize - 1};
        int C_10[4] = {arraySize, 2*arraySize - 1, 0, arraySize - 1};
        int C_11[4] = {arraySize, 2*arraySize - 1, arraySize, 2*arraySize - 1};


        SpecialArray garb0ne = square_matrix_mult_recurse(A, B, A_00, B_00, arraySize);
        SpecialArray garbTwo = square_matrix_mult_recurse(A, B, A_01, B_10, arraySize);
        addMatrices(&garb0ne, &garbTwo, C, C_00, arraySize);
        
        SpecialArray garbThree = square_matrix_mult_recurse(A, B, A_00, B_01, arraySize);
        SpecialArray garbFour = square_matrix_mult_recurse(A, B, A_01, B_11, arraySize);
        addMatrices(&garbThree, &garbFour, C, C_01, arraySize);

        SpecialArray garbFive = square_matrix_mult_recurse(A, B, A_10, B_00, arraySize);
        SpecialArray garbSix = square_matrix_mult_recurse(A, B, A_11, B_10, arraySize);
        addMatrices(&garbFive, &garbSix, C, C_10, arraySize);
        
        SpecialArray garbSeven = square_matrix_mult_recurse(A, B, A_10, B_01, arraySize);
        SpecialArray garbEight = square_matrix_mult_recurse(A, B, A_11, B_11, arraySize);
        addMatrices(&garbSeven, &garbEight, C, C_11, arraySize);
    }
    return C;
}

// 0 is rowStart, 1 is rowEnd, 2 is colStart, 3 is colEnd
SpecialArray square_matrix_mult_strassen(SpecialArray *A, SpecialArray *B, int A_indices[],
                                        int B_indices[], int arraySize){
    //size_t n = A->arrSize;
    SpecialArray C = SpecialArray(arraySize);
    if(arraySize == 1){
        size_t rowA = A_indices[0];
        size_t colA = A_indices[2];
        size_t rowB = B_indices[0];
        size_t colB = B_indices[2];
        C.rows[0][0] = A->rows[rowA][colA]*B->rows[rowB][colB];
    }
    else{
        arraySize /= 2;

        int A_00[4] = {A_indices[0], A_indices[0] + arraySize - 1, A_indices[2],
                        A_indices[2] + arraySize - 1};
        int A_01[4] = {A_indices[0], A_indices[0] + arraySize - 1, A_indices[2] + arraySize, 
                        A_indices[3]};
        int A_10[4] = {A_indices[0] + arraySize, A_indices[1], A_indices[2], 
                        A_indices[2] + arraySize - 1};
        int A_11[4] = {A_indices[0] + arraySize, A_indices[1], A_indices[2] + arraySize,
                        A_indices[3]};
        int B_00[4] = {B_indices[0], B_indices[0] + arraySize - 1, B_indices[2],
                        B_indices[2] + arraySize - 1};
        int B_01[4] = {B_indices[0], B_indices[0] + arraySize - 1, B_indices[2] + arraySize, 
                        B_indices[3]};
        int B_10[4] = {B_indices[0] + arraySize, B_indices[1], B_indices[2], 
                        B_indices[2] + arraySize - 1};
        int B_11[4] = {B_indices[0] + arraySize, B_indices[1], B_indices[2] + arraySize,
                        B_indices[3]};
        
        int C_00[4] = {0, arraySize - 1, 0, arraySize - 1};
        int C_01[4] = {0, arraySize - 1, arraySize, 2*arraySize - 1};
        int C_10[4] = {arraySize, 2*arraySize - 1, 0, arraySize - 1};
        int C_11[4] = {arraySize, 2*arraySize - 1, arraySize, 2*arraySize - 1};

        int defaultArray[4] = {0, arraySize - 1, 0, arraySize -1};

        SpecialArray S1 =  SpecialArray(arraySize);
        subtractMatrices(B, B, S1, B_01, B_11, defaultArray, arraySize);

        SpecialArray S2 =  SpecialArray(arraySize);
        addMatrices(A, A, S2, A_00, A_01, defaultArray, arraySize);

        SpecialArray S3 = SpecialArray(arraySize);
        addMatrices(A, A, S3, A_10, A_11, defaultArray, arraySize);

        SpecialArray S4 =  SpecialArray(arraySize);
        subtractMatrices(B, B, S4, B_10, B_00, defaultArray, arraySize);

        SpecialArray S5 = SpecialArray(arraySize);
        addMatrices(A, A, S5, A_00, A_11, defaultArray, arraySize);

        SpecialArray S6 =  SpecialArray(arraySize);
        addMatrices(B, B, S6, B_00, B_11, defaultArray, arraySize);

        SpecialArray S7 =  SpecialArray(arraySize);
        subtractMatrices(A, A, S7, A_01, A_11, defaultArray, arraySize);

        SpecialArray S8 =  SpecialArray(arraySize);
        addMatrices(B, B, S8, B_10, B_11, defaultArray, arraySize);

        SpecialArray S9 =  SpecialArray(arraySize);
        subtractMatrices(A, A, S9, A_00, A_10, defaultArray, arraySize);

        SpecialArray S10 = SpecialArray(arraySize);
        addMatrices(B, B, S10, B_00, B_01, defaultArray, arraySize);

        SpecialArray P1 = square_matrix_mult_strassen(A, &S1, A_00, defaultArray, arraySize);
        SpecialArray P2 = square_matrix_mult_strassen(&S2, B, defaultArray, B_11, arraySize);
        SpecialArray P3 = square_matrix_mult_strassen(&S3, B, defaultArray, B_00, arraySize);
        SpecialArray P4 = square_matrix_mult_strassen(A, &S4, A_11, defaultArray, arraySize);
        SpecialArray P5 = square_matrix_mult_strassen(&S5, &S6, defaultArray, defaultArray, 
                                                        arraySize);
        SpecialArray P6 = square_matrix_mult_strassen(&S7, &S8, defaultArray, defaultArray, 
                                                        arraySize);
        SpecialArray P7 = square_matrix_mult_strassen(&S9, &S10, defaultArray, defaultArray, 
                                                        arraySize);


        addMatrices(&P5, &P4, C, C_00, arraySize);
        subtractSingleMatrix(&P2, C, defaultArray, C_00, arraySize);
        addSingleMatrix(&P6, C, defaultArray, C_00, arraySize);
        //std::vector<std::vector<float>> C_00 = addMatrices(P5, P4);
        //C_00 = subtractMatrices(C_00, P2);
        //C_00 = addMatrices(C_00, P6);

        addMatrices(&P1, &P2, C, C_01, arraySize);
        //std::vector<std::vector<float>> C_01 = addMatrices(P1, P2);

        addMatrices(&P3, &P4, C, C_10, arraySize);
        //std::vector<std::vector<float>> C_10 = addMatrices(P3, P4);

        addMatrices(&P5, &P1, C, C_11, arraySize);
        subtractSingleMatrix(&P3, C, defaultArray, C_11, arraySize);
        subtractSingleMatrix(&P7, C, defaultArray, C_11, arraySize);
        //std::vector<std::vector<float>> C_11 = addMatrices(P5, P1);
        //C_11 = subtractMatrices(C_11, P3);
        //C_11 = subtractMatrices(C_11, P7);
    }

    return C;
}

SpecialArray generateRandomMatrix(size_t matrixSize){
    
    int nextSize = pow(2, ceil(log(matrixSize)/log(2)));
    SpecialArray* randMatrix = new SpecialArray(nextSize);
    
    for (size_t i = 0; i < matrixSize; i++){
        for (size_t j = 0; j < matrixSize; j++){
            randMatrix->rows[i][j] = rand() % 1000;
        }
    }
    return *randMatrix;
}

void printCSV(std::vector<float> times, std::string fileName){
    std::vector<std::string> labels = {"30", "60", "90", "120", "150", "180", "210", 
    "240", "270", "300", "330", "360", "390"};//, "420", "450", "480", "510", "540", "570", "600"};
    std::ofstream csv;
    csv.open(fileName);
    csv << "input,time\n";
    csv <<"0,0\n";
    for (size_t i = 0; i < labels.size(); i++){
        csv << labels[i] << "," << times[i] << "\n";
    }
    csv.close();
}

void performBruteForceTrials(){
    std::vector<int> nValues = {30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390};

    std::vector<float> times;

    for (size_t i = 0; i < nValues.size(); i++){
        float averateTime = 0.0;
        std::cout << nValues[i] << std::endl <<std::endl;
        for (size_t j = 0; j < TRIALS_PER_N_VALUE; j++){
            SpecialArray A = generateRandomMatrix(nValues[i]);
            SpecialArray B = generateRandomMatrix(nValues[i]);
            int initialA_Indices[] = {0, A.arrSize - 1, 0, A.arrSize- 1};
            int initialB_Indices[] = {0, B.arrSize - 1, 0, B.arrSize - 1};
            auto start = std::chrono::high_resolution_clock::now();
            SpecialArray C = square_matrix_mult(&A, &B, A.arrSize);
            auto end = std::chrono::high_resolution_clock::now();
            //printMatrix(C, C.size());
            //std::cout << "Finished timing" << std::endl <<std::endl;
            averateTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        }
        averateTime /= TRIALS_PER_N_VALUE;
        times.push_back(averateTime);
    }
    
    printCSV(times, "BruteForce.csv");
}

void performStrassenTrials(){
    std::vector<int> nValues = {30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390};

    std::vector<float> times;

    for (size_t i = 0; i < nValues.size(); i++){
        float averateTime = 0.0;
        std::cout << nValues[i] << std::endl <<std::endl;
        for (size_t j = 0; j < TRIALS_PER_N_VALUE; j++){
            SpecialArray A = generateRandomMatrix(nValues[i]);
            SpecialArray B = generateRandomMatrix(nValues[i]);
            int initialA_Indices[] = {0, A.arrSize - 1, 0, A.arrSize- 1};
            int initialB_Indices[] = {0, B.arrSize - 1, 0, B.arrSize - 1};
            auto start = std::chrono::high_resolution_clock::now();
            SpecialArray C = square_matrix_mult_strassen(&A, &B, initialA_Indices, initialB_Indices, A.arrSize);
            auto end = std::chrono::high_resolution_clock::now();
            //std::cout << "Finished timing" << std::endl <<std::endl;
            averateTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        }
        averateTime /= TRIALS_PER_N_VALUE;
        times.push_back(averateTime);
    }
    
    printCSV(times, "Strassen.csv");
}

void performRecurseTrials(){
    std::vector<int> nValues = {30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390};

    std::vector<float> times;

    for (size_t i = 0; i < nValues.size(); i++){
        float averateTime = 0.0;
        std::cout << nValues[i] << std::endl <<std::endl;
        for (size_t j = 0; j < TRIALS_PER_N_VALUE; j++){
            SpecialArray A = generateRandomMatrix(nValues[i]);
            SpecialArray B = generateRandomMatrix(nValues[i]);
            int initialA_Indices[] = {0, A.arrSize - 1, 0, A.arrSize- 1};
            int initialB_Indices[] = {0, B.arrSize - 1, 0, B.arrSize - 1};
            auto start = std::chrono::high_resolution_clock::now();
            SpecialArray C = square_matrix_mult_recurse(&A, &B, initialA_Indices, initialB_Indices, A.arrSize);
            //std::cout << "Finished timing" << std::endl <<std::endl;
            auto end = std::chrono::high_resolution_clock::now();
            averateTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        }
        averateTime /= TRIALS_PER_N_VALUE;
        times.push_back(averateTime);
    }
    
    printCSV(times, "Recurse.csv");
}

int main(){
    performRecurseTrials();
    performStrassenTrials();
    performBruteForceTrials();
    //int size = 4;
    //float *i;
    //i = new float[size*size]{1, 5, 9, 0, 6, 4, 8, 0, 6, 4, 7, 0, 0, 0, 0, 0};
    //float *j;
    //j = new float[size*size]{4, 5, 9, 0, 8, 2, 6, 0, 6, 4, 5, 0, 0, 0, 0, 0};
    
    //SpecialArray *A = new SpecialArray(size, i);
    //SpecialArray *B = new SpecialArray(size, j);
    //float B[size][size] = {{4}};
    //float C[size][size];
    //r |= r >> 1;
    //int nextSize = pow(2, ceil(log(size)/log(2)));

    //std::cout << next;

    //float new_A[nextSize][nextSize];
    //float new_B[nextSize][nextSize];
    //square_matrix_mult(*A, *B, *C, nextSize);
    //int initialA_Indices[] = {0, size - 1, 0, size - 1};
    //int initialB_Indices[] = {0, size - 1, 0, size - 1};
    //SpecialArray C = square_matrix_mult_recurse(A, B, initialA_Indices, initialB_Indices, size);
    //printMatrix(&C, size);
    //SpecialArray D= square_matrix_mult_strassen(A, B, initialA_Indices, initialB_Indices, size);
    //printMatrix(&D, size);
    //A->rows[1][0] = 57;
    //printMatrix(A, size);
    
    /* 44 15 
        56 38 */
    
    /*98 51 84 0 
    
    104 70 118 0 
    
    98 66 113 0 
    
    0 0 0 0 */
    //delete A;
    //delete B;
    return 0;
}