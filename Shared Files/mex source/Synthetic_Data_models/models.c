/*
 * models.c
*/

#include "models.h"
#include "math.h"
#include "stdio.h"


void h_func(double *x, double *lm, double *mu, double *Hl, double *Hn)
{
    double dx = lm[0] - x[0];
    double dy = lm[1] - x[1];
    double d2 = pow(dx,2) + pow(dy,2);
    double d = sqrt(d2);
    
    // Mean
    mu[0] = d;
    mu[1] = pi_to_pi(atan2(dy,dx) - x[2]);
    
    // Jacobian wrt landmark state
    Hl[0] = dx/d;
    Hl[1] = -dy/d2;
    Hl[2] = dy/d;
    Hl[3] = dx/d2;
    
    // Jacobian wrt vehicle state
    Hn[0] = -dx/d;
    Hn[1] = dy/d2;
    Hn[2] = -dy/d;
    Hn[3] = -dx/d2;
    Hn[4] = 0.0;
    Hn[5] = -1.0;
}

void f_func(double *x, double *u, double dt, double *mu, double *Fx, double *Fu)
{
    double s = sin(x[2]);
    double ss = sin(x[2] + u[1]*dt);
    double c = cos(x[2]);
    double cc = cos(x[2] + u[1]*dt);
    double a = u[0]/u[1];
    double a2 = u[0]/pow(u[1],2);
    
    // Mean
    mu[0] = x[0] - a*s + a*ss;
    mu[1] = x[1] + a*c - a*cc;
    mu[2] = x[2] + u[1]*dt;
    
    // Jacobian wrt vehicle state
    Fx[0] = 1.0;
    Fx[1] = 0.0;
    Fx[2] = 0.0;
    Fx[3] = 0.0;
    Fx[4] = 1.0;
    Fx[5] = 0.0;
    Fx[6] = -a*c+a*cc;
    Fx[7] = -a*s+a*ss;
    Fx[8] = 1.0;
    
    // Jacobian wrt control input
    Fu[0] = -1.0/u[1]*s + 1.0/u[1]*ss;
    Fu[1] =  1.0/u[1]*c - 1.0/u[1]*cc;
    Fu[2] = 0.0;
    Fu[3] = a2*s - a2*ss + a*dt*cc;
    Fu[4] = -a2*c + a2*cc + a*dt*ss;
    Fu[5] = dt;
}

void adaptiveDetectionProbability(double *PD, double *xl, double *xn, double FOV_RANGE, double FOV_ANGLE, double DETECTION_PROBABILITY) {
     
    // adaptive detection probability
    double dx = xl[0] - xn[0];
    double dy = xl[1] - xn[1];
    double dx2 = pow(dx,2);
    double dy2 = pow(dy,2);
    
    if (dx2 + dy2 > pow(FOV_RANGE,2))
    {
        *PD = 0.0;
    }
    else if (fabs(pi_to_pi(atan2(dy,dx) - xn[2])) > FOV_ANGLE)
    {
        *PD = 0.0;
    }
    else
    {
        *PD = DETECTION_PROBABILITY;
    }
}


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


void matrixInverse(double *invX, double *det, double *X, int dim)
{   
    double tmp;
    switch(dim){
        case 2:          
            *det = X[0]*X[3] - X[2]*X[1];
            tmp = 1.0 / (*det);
            
            invX[0] =  (tmp) * X[3];
            invX[1] = -(tmp) * X[1];
            invX[2] = -(tmp) * X[2];
            invX[3] =  (tmp) * X[0];
            break;
        case 3:
//             0 3 6 | (0,0)   (0,1)   (0,2)
//             1 4 7 | (1,0)   (1,1)   (1,2)
//             2 5 8 | (2,0)   (2,1)   (2,2)
            
            *det  = X[0] * (X[4]*X[8] - X[7]*X[5]);
            *det -= X[3] * (X[1]*X[8] - X[7]*X[2]);
            *det += X[6] * (X[1]*X[5] - X[4]*X[2]);
            
            tmp = 1.0 / (*det);
            
            invX[0] = (tmp) * (X[4] * X[8] - X[5] * X[7]);
            invX[3] = (tmp) * (X[6] * X[5] - X[3] * X[8]);
            invX[6] = (tmp) * (X[3] * X[7] - X[6] * X[4]);
            invX[1] = (tmp) * (X[7] * X[2] - X[1] * X[8]);
            invX[4] = (tmp) * (X[0] * X[8] - X[6] * X[2]);
            invX[7] = (tmp) * (X[1] * X[6] - X[0] * X[7]);
            invX[2] = (tmp) * (X[1] * X[5] - X[2] * X[4]);
            invX[5] = (tmp) * (X[2] * X[3] - X[0] * X[5]);
            invX[8] = (tmp) * (X[0] * X[4] - X[1] * X[3]);
            break;
    }
}

void choleskyDecomposition(double *A, double *B)
{
    B[0] = sqrt(A[0]);
    B[1] = A[1]/B[0];
    B[2] = A[2]/B[0];
    
    B[3] = 0.0;
    B[4] = sqrt(A[4] - pow(B[1],2));
    B[5] = (A[5] - B[2]*B[1])/B[4];
    
    B[6] = 0.0;
    B[7] = 0.0;
    B[8] = sqrt(A[8] - pow(B[2],2) - pow(B[5],2));
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

