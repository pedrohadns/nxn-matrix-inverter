#include <iostream>
#include <stdio.h>

float **allocate(int n){ //Allocates the necessary space for a certain nxn matrix 'u' memory
    float **u = new float*[n]; //Points to each row of the matrix
    for(int i = 0; i < n; i++){ 
        u[i] = new float[n]; //Each ith row has a pointer of size n to each element u_i, j of the matrix
    }

    return u;
}

void deallocate(float **u, int n){ //Deallocates a certain nxn matrix 'u' from memory
    for(int i = 0; i < n; i++){
        delete[] u[i];
    }
    delete[] u;
}

void printMatrix(float **u, int n, const char *s1 = "", const char *s2 = ""){ //Prints a nxn matrix 'u' in a formatted way, s1 and s2 are the elements that surround the matrix, e.g (1, 2 ,3) '(' and ')' are s1 and s2, respectively
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(j == 0){
                printf("%s%.3f ", s1, u[i][j]);
            }
            else if(j == n - 1){
                printf("%.3f%s\n", u[i][j], s2);
            }
            else{
                printf("%.3f ", u[i][j]);
            }
        }
    }
    printf("\n");
}

float **input(int n){ //Asks the user to input each element of a nxn matrix
    float **m = allocate(n);

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("Type the element m_%d, %d of the matrix: ", i + 1, j + 1);
            scanf("%f", &m[i][j]);
        }
    }

    return m;
}

bool isThereZero(float **u, int n){ //Verifies if a matrix has 0 as its interior element(s)
    if(n == 1){
        return false; //1x1 matrixes don't have to be verified, as they are already the final result of the condensation
    }
        
    else{
        bool anyZeros = false;
        if(n == 2){
            return anyZeros;
        }
        else{
            for(int i = 1; i < n - 1; i++){
                for(int j = 1; j < n - 1; j++){
                    anyZeros = u[i][j] == 0;
                    if(anyZeros == true){break;}
                }
                if(anyZeros == true){break;}
            }
            return anyZeros;
        }
    }
}

void addRows(float **u, int n, int sourceRow, int beingAddedRow){ //Adds two rows of a matrix
    for(int j = 0; j < n; j++){
        u[sourceRow][j] += u[beingAddedRow][j];
    }
}

float **eliminateZero(float **u, int n){ //Eliminates zeros in interior by adding rows of a matrix in order to avoid division by 0 during Dodgson's condensation
    float **non_Zero = allocate(n);

    for(int i = 0; i < n; i++){ //Creates a matrix with all elements of the original that will have its zeros removed
        for(int j = 0; j < n; j++){
            non_Zero[i][j] = u[i][j];
        }
    }
    int center;
    if(n % 2 == 0){
        center = (n / 2) - 1;
    }
    else{
        center = (n - 1) / 2;
    }
    
    while(isThereZero(non_Zero, n)){
        for(int i = 1; i < n - 1; i++){
            for(int j = 1; j < n - 1; j++){
                for(int count = 1; count < n - 1; count++){
                    if(non_Zero[i][j] == 0){
                        if(non_Zero[i + count][j] != 0){
                            addRows(non_Zero, n, i, i + count);
                            break;
                        }
                        else if((non_Zero[i - count][j] != 0) && (count <= center)){
                            addRows(non_Zero, n, i, i - count);
                            break;
                        }
                    }
                }
            }
        }
    }

    return non_Zero;
}


float** determinant(float **A, float **B, float **C, int n, int step = 0){ //Calculates the determinant of a matrix using Dodgson's condensation method
    A = eliminateZero(A, n);
    
    B = allocate(n - 1);
    if(n != 2 && n != 1){
        for(int i = 0; i < (n - 1); i++){
            for(int j = 0; j < (n - 1); j++){
                float Temporary[2][2] = {{A[i][j],     A[i][j + 1]}, 
                                         {A[i + 1][j], A[i + 1][j + 1]}};
                float detTemporary = ((Temporary[0][0] * Temporary[1][1]) - (Temporary[0][1] * Temporary[1][0]));

                if(step != 0){
                    B[i][j] = detTemporary / C[i + 1][j + 1];
                }
                else{
                    B[i][j] = detTemporary;
                }
            }
        }

        if((step % 2) == 0){
            C = allocate(n);
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    C[i][j] = A[i][j];
                }
            }
        }
        step++;
        
        determinant(B, B, A, n - 1, step);
        
    }

    else if(n == 2){
        if(step != 0){
            float **detA = allocate(1);
            detA[0][0] = ((A[0][0] * A[1][1]) - (A[0][1] * A[1][0]));
            step++;

            determinant(detA, B, C, 1, step);
        }

        else{
            float **detA = allocate(1);
            detA[0][0] = ((A[0][0] * A[1][1]) - (A[0][1] * A[1][0]));

            return detA;
        }
    }

    else if(n == 1){
        if(step != 0){
            float **detA = allocate(1);
        
            detA[0][0] = A[0][0] / C[1][1];

            return detA;
        }
        else{
            return A;
        }
    }
}

float **minorMatrix(float **m, int n){ //Calculates the minor matrix
    float **T = allocate(n - 1); //Temporary matrix that receives the elements from the rows that aren't disregarded
    float **X = allocate(n);
    float **B = allocate(n); float **C = allocate(n);

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            B[i][j] = 1; C[i][j] = 1;
        }
    }

    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            int sub_i = 0;
            for (int i = 0; i < n; i++) {
                if (i == row) continue;
                int sub_j = 0;
                for (int j = 0; j < n; j++) {
                    if (j == col) continue;
                    T[sub_i][sub_j] = m[i][j];
                    sub_j++;
                }
                sub_i++;
            }

            X[row][col] = determinant(T, B, C, n - 1)[0][0];
             
        }
    }

    deallocate(T, n - 1); deallocate(B, n); deallocate(C, n);

    return X;
}

float **matrixCofactors(float **m, int n){ //Calculates the matrix of cofactors
    float **cofactor = allocate(n);

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(((i + j) % 2) == 0){
                cofactor[i][j] = m[i][j];
            }
            else{
                cofactor[i][j] = -m[i][j];
            }
        }
    }

    return cofactor;
}

float **matrixAdjoint(float **m, int n){ //Calculates an adjoint matrix
    float **adjoint = allocate(n);

        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                adjoint[i][j] = m[j][i];
            }
        }

    return adjoint;
}

float **multiplyScalar(float **m, int n, float k){ //Multiplies a nxn matrix by a scalar
    float **multiplied = allocate(n);

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            multiplied[i][j] = m[i][j] * k;
        }
    }

    return multiplied;
}

int main(){
    int n;
    printf("Type the matrix size: ");
    scanf("%d", &n);

    //int n = 6;

    float **B = allocate(n); float **C = allocate(n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            B[i][j] = 1; C[i][j] = 1;
        }
    }

    float **m = input(n);

    /*float **m = allocate(n);
    m[0][0] = 6; m[0][1] = 1; m[0][2] = 9; m[0][3] = 1; m[0][4] = 0; m[0][5] = 3;
    m[1][0] = 9; m[1][1] = 9; m[1][2] = 1; m[1][3] = 9; m[1][4] = 4; m[1][5] = 9;
    m[2][0] = 7; m[2][1] = 4; m[2][2] = 4; m[2][3] = 3; m[2][4] = 9; m[2][5] = 4;
    m[3][0] = 7; m[3][1] = 3; m[3][2] = 9; m[3][3] = 8; m[3][4] = 4; m[3][5] = 0;
    m[4][0] = 8; m[4][1] = 6; m[4][2] = 0; m[4][3] = 0; m[4][4] = 0; m[4][5] = 4;
    m[5][0] = 5; m[5][1] = 5; m[5][2] = 5; m[5][3] = 5; m[5][4] = 7; m[5][5] = 8;*/

    printf("\n---------------MATRIX---------------n");
    printMatrix(m, n);

    float detM = determinant(m, B, C, n)[0][0];

    if(detM != 0){
        printf("\n---------------MINOR---------------n");
        float **X = minorMatrix(m, n);
        printMatrix(X, n, "[", "]");

        printf("\n---------------COFACTOR---------------n");
        float **cofactor = matrixCofactors(X, n);
        printMatrix(cofactor, n, "[", "]");
        deallocate(X, n);

        printf("\n---------------ADJOINT---------------n");
        float **adjoint = matrixAdjoint(cofactor, n);
        printMatrix(adjoint, n, "[", "]");
        deallocate(cofactor, n);

        printf("\nDeterminant: %f\n", detM);

        printf("\n---------------INVERSE---------------\n");
        float **multiplied = multiplyScalar(adjoint, n, 1 / detM);
        deallocate(adjoint, n);
        printMatrix(multiplied, n, "(", ")");

        deallocate(m, n); deallocate(B, n); deallocate(C, n); deallocate(multiplied, n);
        
        return 0;
    }
    else{
        printf("Matrix can't be inverted: determinant = 0.\n");
        deallocate(m, n); deallocate(B, n); deallocate(C, n);

        return 1;
    }
}