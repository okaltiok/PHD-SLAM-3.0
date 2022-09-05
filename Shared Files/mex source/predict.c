//     This function implements the PHD prediction step including the birth 
//     process that creates new landmarks the propagation of the vehicle
//     according to the kinematic model and control inputs
// 
//     Input:
//        obj        - struct that represent particle n of the PHD-SLAM density at time k-1
//        u          - {2 x T} control matrix
//        t          - {1 x T} time vector
//        params     - simulation parameters
//     
//     Output:
//        obj        - struct that represent particle n of the PHD-SLAM density at time k | k-1
//     
//     Author   : Ossi Kaltiokallio
//                Tampere University, Department of Electronics and
//                Communications Engineering
//                Korkeakoulunkatu 1, 33720 Tampere
//                ossi.kaltiokallio@tuni.fi
//     Last Rev : 26/8/2022
//     Tested   : '9.8.0.1359463 (R2020a) Update 1'
//     
//     Copyright notice: You are free to modify, extend and distribute 
//        this code granted that the author of the original code is 
//        mentioned as the original author of the code.


#include "mex.h"
#include "matrix.h"
#include "string.h"
#include "math.h"
#include "models.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*
     * nrhs     Number of input (right-side) arguments, or the size of the prhs array.
     * prhs     Array of input arguments.
     */
     
    /* duplicate input to output */
    plhs[0] = mxDuplicateArray(prhs[0]);
    
    /* variable declarations of field pointer */
    mxArray* birth_ptr = mxGetField(plhs[0], 0, "birth_y"); 
    
    /* variable declarations of output variables */
    double* xn = mxGetDoubles(mxGetField(plhs[0], 0, "xn"));
    double* Pn = mxGetDoubles(mxGetField(plhs[0], 0, "Pn"));
    double* birth_y = mxGetDoubles(birth_ptr);
    double* u = mxGetDoubles(prhs[1]);      // 5 x ncols measurement matrix
    double* t = mxGetDoubles(prhs[2]);      // 5 x ncols measurement matrix
    
    /* variable declarations of constants */
    double* Qu = mxGetDoubles(mxGetField(prhs[3], 0, "Qu"));
    double* Qn = mxGetDoubles(mxGetField(prhs[3], 0, "Qn"));
    double* R = mxGetDoubles(mxGetField(prhs[3], 0, "R"));
    double P_B = mxGetScalar(mxGetField(prhs[3], 0, "P_B"));
    double threshold = mxGetScalar(mxGetField(prhs[3], 0, "eta_threshold"));
    
    /* dimension declarations */
    int xn_dim = mxGetScalar(mxGetField(prhs[3], 0, "xn_dim"));
    int xl_dim = mxGetScalar(mxGetField(prhs[3], 0, "xl_dim"));
    int u_dim = mxGetScalar(mxGetField(prhs[3], 0, "u_dim"));
    int h_dim = mxGetScalar(mxGetField(prhs[3], 0, "h_dim"));
    int m_k = mxGetN(birth_ptr);            // number of measurements
    int n_k = mxGetN(mxGetField(prhs[0], 0, "xl"));      // number of landmarks

    /* create new landmarks */
    if (m_k > 0)
    {
        /* variable declarations of output variables */
        mxArray* xl_ptr = mxGetField(plhs[0], 0, "xl");
        mxArray* Pl_ptr = mxGetField(plhs[0], 0, "Pl");
        mxArray* eta_ptr = mxGetField(plhs[0], 0, "eta");
        mxArray* eta_threshold_ptr = mxGetField(plhs[0], 0, "eta_threshold");
        
        double* xl = mxGetDoubles(xl_ptr);
        double* Pl = mxGetDoubles(Pl_ptr);
        double* eta = mxGetDoubles(eta_ptr);
        double* eta_threshold = mxGetDoubles(eta_threshold_ptr);
        
        /* define dimensions of landmark parameters */
        mwSize xl_dims[2] = {xl_dim,n_k + m_k};
        mwSize Pl_dims[3] = {xl_dim,xl_dim,n_k + m_k};
        mwSize eta_dims[2] = {1,n_k + m_k};
        
        if (n_k == 0)
        {
            xl_dims[0] = xl_dim;
            Pl_dims[0] = xl_dim;
            Pl_dims[1] = xl_dim;
            eta_dims[0] = 1;
        }

        /* reallocate memory for new landmarks */
        mxSetDimensions(xl_ptr, xl_dims, 2);
        xl = mxRealloc(xl,(n_k+m_k)*xl_dim*sizeof(double));
        mxSetPr(xl_ptr,xl);
        
        mxSetDimensions(eta_ptr, eta_dims, 2);
        eta = mxRealloc(eta,(n_k+m_k)*sizeof(double));
        mxSetPr(eta_ptr,eta);
        
        mxSetDimensions(eta_threshold_ptr, eta_dims, 2);
        eta_threshold = mxRealloc(eta_threshold,(n_k+m_k)*sizeof(double));
        mxSetPr(eta_threshold_ptr,eta_threshold);

        mxSetDimensions(Pl_ptr, Pl_dims, 3);
        Pl = mxRealloc(Pl,(n_k+m_k)*xl_dim*xl_dim*sizeof(double));
        mxSetPr(Pl_ptr,Pl);
    
        xl += xl_dim*n_k;
        Pl += xl_dim*xl_dim*n_k;
        eta += n_k;
        eta_threshold += n_k;

        double G[4], GR[4];
        double r, b, s, c;
        
        /* create new landmark for every measurement */
        for (int i = 0; i < m_k; i++)
        {
            r = *(birth_y++);
            b = *(birth_y++);
            s = sin(xn[2]+b);
            c = cos(xn[2]+b);
            
            // mean
            xl[0] = xn[0] + r*c;
            xl[1] = xn[1] + r*s;
            
            // Jacobian
            G[0] = c;
            G[1] = s;
            G[2] = -r*s;
            G[3] = r*c;
            
            /* compute covariance Pl = G*R*G' */
            matrixMultiply(xl_dim, xl_dim, xl_dim, xl_dim, 0, 0, G, R, GR);
            matrixMultiply(xl_dim, xl_dim, xl_dim, xl_dim, 1, 0, GR, G, Pl);
            
            *eta = log(P_B);
            *eta_threshold = threshold;
            
            xl += xl_dim;
            Pl += xl_dim*xl_dim;
            eta++;
            eta_threshold++;
        }
        // deallocate memory
        birth_y -= m_k*h_dim;
        const mwSize empty_ndim = 2;
        const mwSize empty_dims[2] = {0, 0};
        mxSetDimensions(birth_ptr, empty_dims, empty_ndim);
        birth_y = mxRealloc(birth_y,0);
        mxSetPr(birth_ptr,birth_y);
    }

    /* propagate vehicle */
    int T = mxGetN(prhs[2]);      
    if (T > 1)
    {
        for (int i=0; i < xn_dim*xn_dim; i++)
            Pn[i] += Qn[i];

        double dt, mu[3], Fx[9], FXP[9], Fu[6], FUQ[6];
        for (int i = 1; i < T; i++)
        {
            dt = t[i] - t[i-1];
            
            /* compute mean and Jacobians */
            f_func(xn, u, dt, mu, Fx, Fu);

            /* compute Pn = Fx*Pn*Fx' + Pn*/
            matrixMultiply(xn_dim, xn_dim, xn_dim, xn_dim, 0, 0, Fx, Pn, FXP);
            matrixMultiply(xn_dim, xn_dim, xn_dim, xn_dim, 1, 0, FXP, Fx, Pn);

            /* compute Pn = Fu*Q*Fu' + Pn*/
            matrixMultiply(xn_dim, xl_dim, xl_dim, xl_dim, 0, 0, Fu, Qu, FUQ);
            matrixMultiply(xn_dim, xl_dim, xn_dim, xl_dim, 1, 1, FUQ, Fu, Pn);

            memcpy(xn, mu, xn_dim*sizeof(double) );
            u += u_dim;    
        }
    }
}




