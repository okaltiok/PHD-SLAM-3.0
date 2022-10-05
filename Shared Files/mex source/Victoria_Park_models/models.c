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

void f_func(double *m, double *u, double dt, double *mu, double *Fx, double *Fu)
{
    double L = 2.83;
    double H = 0.75; 
    double a = 3.78;
    double b = 0.5;
    
    // compute help variables
    double v_c = u[0]/(1-H/L*tan(u[1]));
    double sc = a*sin(m[2]) + b*cos(m[2]);
    double cs = a*cos(m[2]) - b*sin(m[2]);
    double tmp1 = cos(m[2]) - sc*tan(u[1])/L;
    double tmp2 = sin(m[2]) + cs*tan(u[1])/L;
    double tmp3 = dt*v_c/(L*pow(cos(u[1]),2));

    // Mean
    mu[0] = m[0] + dt*v_c*tmp1;
    mu[1] = m[1] + dt*v_c*tmp2;
    mu[2] = m[2] + dt*v_c*tan(u[1])/L;
    
    // Jacobian wrt vehicle state
    Fx[0] = 1.0;
    Fx[1] = 0.0;
    Fx[2] = 0.0;
    Fx[3] = 0.0;
    Fx[4] = 1.0;
    Fx[5] = 0.0;
    Fx[6] = -dt*v_c*tmp2;
    Fx[7] = dt*v_c*tmp1;
    Fx[8] = 1.0;
    
    // Jacobian wrt control input
    Fu[0] = dt*tmp1;
    Fu[1] = dt*tmp2;
    Fu[2] = dt*tan(u[1])/L;
    Fu[3] = -sc*tmp3;
    Fu[4] = cs*tmp3;
    Fu[5] = tmp3;
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
        *PD = (1 - sqrt(dx2 + dy2)/FOV_RANGE)*DETECTION_PROBABILITY;
    }
}