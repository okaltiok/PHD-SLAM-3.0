//     This function computes the vehicle and map estimates
// 
//     Input:
//        obj     - struct that representthe PHD-SLAM density
//        params  - simulation parameters
//     
//     Output:
//        est      - struct that contains the estimates
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
#include "math.h"
#include <string.h>


#define NUMBER_OF_FIELDS (sizeof(field_names) / sizeof(*field_names))


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    
    /* Create output struct */
    const char* field_names[] = {"m_hat", "P_hat", "mu_hat", "C_hat", "active"};
    mwSize struct_dims[2] = {1, 1};
    plhs[0] = mxCreateStructArray(2, struct_dims, NUMBER_OF_FIELDS, field_names);
    
    /* dimension declarations */
    int xn_dim = mxGetScalar(mxGetField(prhs[1], 0, "xn_dim"));
    int xl_dim = mxGetScalar(mxGetField(prhs[1], 0, "xl_dim"));
 
    /* get number of particles and find index of maximum weight*/
    mwSize N = mxGetNumberOfElements(prhs[0]);   
    mwIndex jstruct, jmax;
    double w, w_max;
    for(jstruct = 0; jstruct < N; jstruct++) {
        w = *mxGetDoubles(mxGetFieldByNumber(prhs[0], jstruct, mxGetFieldNumber(prhs[0], "w")));
        if(jstruct == 0 || w > w_max)
        {
            w_max = w;
            jmax = jstruct;
        }
    }

    /* get fields corresponding to maximum weight*/
    double* xn = mxGetDoubles(mxGetField(prhs[0], jmax, "xn"));
    double* Pn = mxGetDoubles(mxGetField(prhs[0], jmax, "Pn"));
    double* xl = mxGetDoubles(mxGetField(prhs[0], jmax, "xl"));
    double* Pl = mxGetDoubles(mxGetField(prhs[0], jmax, "Pl"));
    double* eta = mxGetDoubles(mxGetField(prhs[0], jmax, "eta"));
    double* eta_threshold = mxGetDoubles(mxGetField(prhs[0], jmax, "eta_threshold"));
    double* xl_p = mxGetDoubles(mxGetField(prhs[0], jmax, "xl_p"));
    double* Pl_p = mxGetDoubles(mxGetField(prhs[0], jmax, "Pl_p"));
    double* eta_p = mxGetDoubles(mxGetField(prhs[0], jmax, "eta_p"));
    double* eta_threshold_p = mxGetDoubles(mxGetField(prhs[0], jmax, "eta_threshold_p"));
    int xl_N = mxGetN(mxGetField(prhs[0], jmax, "eta"));
    int xlp_N = mxGetN(mxGetField(prhs[0], jmax, "eta_p"));
    
    /* compute number of estimates within the FoV */
    int nk = 0, nk_p = 0; 
    for(int k = 0; k < xl_N; k++) {
        if(eta[k] > log(eta_threshold[k]))
            nk++; 
    }
    
    /* compute number of estimates outside the FoV */
    for(int k = 0; k < xlp_N; k++) {
        if(eta_p[k] > log(eta_threshold_p[k]))
            nk_p++;
    }
    
    /* define output variables */
    const mwSize ndim = 3;
    const mwSize dims[3] = {xl_dim,xl_dim,nk + nk_p};
    
    mxArray* m_ptr = mxCreateDoubleMatrix(xn_dim, 1, mxREAL);
    mxArray* P_ptr = mxCreateDoubleMatrix(xn_dim, xn_dim, mxREAL);
    mxArray* mu_ptr = mxCreateDoubleMatrix(xl_dim, nk + nk_p, mxREAL);
    mxArray* C_ptr = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    mxArray* active_ptr = mxCreateDoubleMatrix(1, nk + nk_p, mxREAL);
    
    double* m_hat = mxGetDoubles(m_ptr);
    double* P_hat = mxGetDoubles(P_ptr);
    double* mu_hat = mxGetDoubles(mu_ptr);
    double* C_hat = mxGetDoubles(C_ptr);
    double* active = mxGetDoubles(active_ptr);
    
    /* copy data */
    memcpy(m_hat, xn, xn_dim*sizeof(double) );
    memcpy(P_hat, Pn, xn_dim*xn_dim*sizeof(double) );
    
    for(int k = 0; k < xl_N; k++) {
        if(eta[k] > log(eta_threshold[k]))
        {
            memcpy(mu_hat, xl + k*xl_dim, xl_dim*sizeof(double) );
            memcpy(C_hat, Pl + k*xl_dim*xl_dim, xl_dim*xl_dim*sizeof(double) );
            *active = 1;
            mu_hat += xl_dim;
            C_hat += xl_dim*xl_dim;
            active++;
        }
    }
    
    for(int k = 0; k < xlp_N; k++) {
        if(eta_p[k] > log(eta_threshold_p[k]))
        {
            memcpy(mu_hat, xl_p + k*xl_dim, xl_dim*sizeof(double) );
            memcpy(C_hat, Pl_p + k*xl_dim*xl_dim, xl_dim*xl_dim*sizeof(double) );
            *active = 0;
            mu_hat += xl_dim;
            C_hat += xl_dim*xl_dim;
            active++;
        }
    }

    /* set fields of output struct */
    mxSetField(plhs[0], 0, field_names[0], m_ptr);
    mxSetField(plhs[0], 0, field_names[1], P_ptr);
    mxSetField(plhs[0], 0, field_names[2], mu_ptr);
    mxSetField(plhs[0], 0, field_names[3], C_ptr);
    mxSetField(plhs[0], 0, field_names[4], active_ptr);
}
