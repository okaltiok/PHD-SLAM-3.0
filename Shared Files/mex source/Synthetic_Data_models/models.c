/*
 * models.c
*/

#include "linearAlgebra.h"
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