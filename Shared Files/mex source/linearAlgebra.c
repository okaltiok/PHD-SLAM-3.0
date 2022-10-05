/*
 * linearAlgebra.c
*/

#include "linearAlgebra.h"
#include "math.h"

double pi_to_pi(double angle) {
    angle = fmod(angle, 2*M_PI);
    if (angle > M_PI)
    {
        angle -= 2*M_PI;
    }
    else if (angle < -M_PI)
    {
        angle += 2*M_PI;
    }
    return angle;
}

bool choleskyFactorization(double *L, double *A, int dim)
{   
    double tmp;
    for (int i = 0; i < dim*dim; i++) 
        L[i] = 0.0;
    
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            for (int k = 0; k < j; k++)
                sum += L[i+k*dim] * L[j+k*dim];

            if (i == j)
            {
                tmp = A[i+i*dim] - sum;
                if (tmp <= 0)
                    return false;
                
                L[i+j*dim] = sqrt(tmp);
            }
            else
            {
                L[i+j*dim] = (A[i+j*dim] - sum) / L[j+j*dim] ;   
            }
        }
    }
    return true;
}

void matrixMultiply(int rA, int cA, int rB, int cB, int transpose, int opt, double *A, double *B, double *C)
{
    if (transpose == 1)
    {
        for (int i = 0; i < rA; i++) {
            for (int j = 0; j < rB; j++) {
                double sum = 0.0;
                for (int k = 0; k < cA; k++)
                {
                    sum = sum + A[i + rA * k] * B[j + rB * k];   
                }
                switch(opt){
                    case 0:
                        C[i + rA * j] = sum;
                        break;
                    case 1:
                        C[i + rA * j] += sum;
                        break;
                    case -1:
                        C[i + rA * j] -= sum;
                        break;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < rA; i++) {
            for (int j = 0; j < cB; j++) {
                double sum = 0.0;
                for (int k = 0; k < cA; k++)
                {
                    sum = sum + A[i + rA * k] * B[j * rB + k];
                }
                switch(opt){
                    case 0:
                        C[i + rA * j] = sum;
                        break;
                    case 1:
                        C[i + rA * j] += sum;
                        break;
                    case -1:
                        C[i + rA * j] -= sum;
                        break;
                }
            }
        }
    }  
}

bool matrixInverse(double *invX, double *det, double *X, int dim)
{   
    switch(dim){
        case 1:    
        {
            *det = X[0];
            
            if (*det == 0)
                return false;
                        
            invX[0]  = 1.0 / (*det);
            return true;
        }
        case 2:    
        {
            *det = X[0]*X[3] - X[2]*X[1];
            
            if (*det == 0)
                return false;
            
            double tmp = 1.0 / (*det);
            
            invX[0] =  (tmp) * X[3];
            invX[1] = -(tmp) * X[1];
            invX[2] = -(tmp) * X[2];
            invX[3] =  (tmp) * X[0];
            return true;
        }
        case 3:
        {
            *det =  X[0] * (X[4]*X[8] - X[7]*X[5]);
            *det -= X[3] * (X[1]*X[8] - X[7]*X[2]);
            *det += X[6] * (X[1]*X[5] - X[4]*X[2]);
            
            if (*det == 0)
                return false;
            
            double tmp = 1.0 / (*det);
            
            invX[0] = (tmp) * (X[4] * X[8] - X[5] * X[7]);
            invX[3] = (tmp) * (X[6] * X[5] - X[3] * X[8]);
            invX[6] = (tmp) * (X[3] * X[7] - X[6] * X[4]);
            invX[1] = (tmp) * (X[7] * X[2] - X[1] * X[8]);
            invX[4] = (tmp) * (X[0] * X[8] - X[6] * X[2]);
            invX[7] = (tmp) * (X[1] * X[6] - X[0] * X[7]);
            invX[2] = (tmp) * (X[1] * X[5] - X[2] * X[4]);
            invX[5] = (tmp) * (X[2] * X[3] - X[0] * X[5]);
            invX[8] = (tmp) * (X[0] * X[4] - X[1] * X[3]);
            return true;
        }
    }
}
