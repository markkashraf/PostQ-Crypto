#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TK_K 3
#define TK_N 4
#define TK_Q 97



static void naivePolynomialMultiply(short *res, const short *poly1, const short *poly2, int acc) {
    short intmd[TK_N];

    intmd[0] = (poly1[0]*poly2[0] + poly1[3]*poly2[1] + poly1[2]*poly2[2] + poly1[1]*poly2[3]) % TK_Q;
    intmd[1] = (poly1[1]*poly2[0] + poly1[0]*poly2[1] + poly1[3]*poly2[2] + poly1[2]*poly2[3]) % TK_Q;
    intmd[2] = (poly1[2]*poly2[0] + poly1[1]*poly2[1] + poly1[0]*poly2[2] + poly1[3]*poly2[3]) % TK_Q;
    intmd[3] = (poly1[3]*poly2[0] + poly1[2]*poly2[1] + poly1[1]*poly2[2] + poly1[0]*poly2[3]) % TK_Q;

    for (int i = 0; i < TK_N; ++i) {
        res[i] = acc ? (res[i] + intmd[i]) % TK_Q : intmd[i];
    }
}



void addPolynomials(short *result, short *poly1, short *poly2){
    for(int index = 0; index < TK_N; index++){
        result[index] = (poly1[index] + poly2[index]) % TK_Q;
    }
}



// p = a.b
void matrixVectorMultiply(short *p, short *a, short *b){
    short r_conv[TK_K*TK_K*TK_N];
    // convolute vertor on matrix
    for (int i=0; i < TK_K; i++) {
            naivePolynomialMultiply(&r_conv[(i*TK_K*TK_N) + (0*TK_N)], &a[(i*TK_K*TK_N) + (0*TK_N)], &b[0*TK_N], 0);
            naivePolynomialMultiply(&r_conv[(i*TK_K*TK_N) + (1*TK_N)], &a[(i*TK_K*TK_N) + (1*TK_N)], &b[1*TK_N], 0);
            naivePolynomialMultiply(&r_conv[(i*TK_K*TK_N) + (2*TK_N)], &a[(i*TK_K*TK_N) + (2*TK_N)], &b[2*TK_N], 0);
    }
    // sum r_conv on axis 0
        addPolynomials(&p[0*TK_N], &r_conv[0*TK_N], &r_conv[1*TK_N]);
        addPolynomials(&p[0*TK_N], &p[0*TK_N], &r_conv[2*TK_N]);

        addPolynomials(&p[1*TK_N], &r_conv[3*TK_N], &r_conv[4*TK_N]);
        addPolynomials(&p[1*TK_N], &p[1*TK_N], &r_conv[5*TK_N]);

        addPolynomials(&p[2*TK_N], &r_conv[6*TK_N], &r_conv[7*TK_N]);
        addPolynomials(&p[2*TK_N], &p[2*TK_N], &r_conv[8*TK_N]);
}


void generateToyPolynomials(short *matrix, short *t, short *s) {
    short rand_noise[TK_K * TK_N];
    srand((unsigned int)time(NULL));

    // Fill A with uniformly random numbers mod q
    for (int i = 0; i < TK_K * TK_K * TK_N; ++i) {
        matrix[i] = rand() % TK_Q;
    }

    // Fill s and e with small random numbers
    for (int i = 0; i < TK_K * TK_N; ++i) {
        s[i] = (rand() & 3)%TK_Q;
        rand_noise[i] = (rand() & 3)%TK_Q;
    }

    // t = A.s + e
    matrixVectorMultiply(t, matrix, s);  // Compute A.s
    addPolynomials(t, t, rand_noise); // Add e

}

void swapRows(short *matrix, int row1, int row2) {
    int startIndex1 = row1 * TK_N;
    int startIndex2 = row2 * TK_N;

    for (int i = 0; i < TK_N; i++) {
        short temp = matrix[startIndex1 + i];
        matrix[startIndex1 + i] = matrix[startIndex2 + i];
        matrix[startIndex2 + i] = temp;
    }
}

void transposeMatrix(short *matrix) {
    swapRows(matrix, 1, 3);
    swapRows(matrix, 5, 7);
    swapRows(matrix, 2, 6);
}



void toyEncrypt( short *A,  short *t, int plain, short *u, short *v){

    short vecR[TK_K * TK_N];
    short vecE1[TK_K * TK_N];
    short vecE2[TK_N];
    short msgBits[TK_N];
    // Fill r, e1, and e2 with small random numbers
    for (int i = 0; i < TK_K * TK_N; ++i) {
        vecR[i] = (rand() & 3)%TK_Q;
        vecE1[i] = (rand() & 3)%TK_Q;
    }
    for (int i = 0; i < TK_N; ++i) {
        vecE2[i] =(rand() & 3)%TK_Q;
        int val = rand() & 3; // Uniform distribution
        vecE2[i] = (val & 1) - ((val >> 1) & 1); // Binomial distribution
        vecE2[i] %= TK_Q;
    }
    // u = A_transpose . r + e1
    transposeMatrix(A);
    matrixVectorMultiply(u, A,vecR);
    addPolynomials(u,u,vecE1);
    for (int i = 0; i < TK_K; ++i) {
        naivePolynomialMultiply(msgBits, &t[i * TK_N], &vecR[i * TK_N], i != 0);
    }
     addPolynomials(v,msgBits,vecE2);
    for (int i = 0; i < TK_N; ++i) {
        v[i] = (v[i] + ((plain >> i) & 1) * (TK_Q / 2))% TK_Q;
    }
}


int toyDecrypt(const short *secretS, const short *vectorU, const short *vectorV) {
    short polyP[TK_N];
    int plaintext = 0;

    // Compute p = v - s . u
    for (int i = 0; i < TK_K; ++i) {
        naivePolynomialMultiply(polyP, &secretS[i * TK_N], &vectorU[i * TK_N], i != 0);
    }
    for (int i = 0; i < TK_N; ++i) {
        polyP[i] = (vectorV[i] - polyP[i] + TK_Q) % TK_Q; // Adjust for negative values
    }

    // Recover plaintext
    for (int i = 0; i < 4; ++i) {
        short x = polyP[i];
        // Adjust if val is > TK_Q/2
        if (x > TK_Q / 2) {
            x -= TK_Q;
        }
        // Determine the bit value based on the magnitude of val
        int bit = abs(x) > TK_Q / 4;
        plaintext |= bit << i;
    }

    return plaintext;
}



int main(){
    short A[TK_K*TK_K*TK_N];
    short t[TK_K*TK_N];
    short s[TK_K*TK_N];
    short u[TK_K*TK_N];
    short v[TK_N];
    generateToyPolynomials(A, t, s);
    toyEncrypt(A, t, 15, u, v);
    int plain;
    plain = toyDecrypt(s, u, v);

    printf("%d", plain);
}