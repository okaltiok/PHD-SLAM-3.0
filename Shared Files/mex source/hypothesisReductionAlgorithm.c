//     This function performs GM component pruning and merging, and
//     determines the components that are active and passive. An active
//     component is within the FOV or close to it, and it is passive 
//     otherwise.
// 
//     Input:
//        obj        - a struct that represent the updated particle (i) of the PHD-SLAM density at time k
//        params     - simulation parameters
// 
//     Output:
//        obj        - obj with pruned and merged GM components
//     
//     Author   : Ossi Kaltiokallio
//                Tampere University, Department of Electronics and
//                Communications Engineering
//                Korkeakoulunkatu 1, 33720 Tampere
//                ossi.kaltiokallio@tuni.fi
//     Last Rev : 30/3/2022
//     Tested   : Matlab version 9.8.0.1359463 (R2020a) Update 1
//     
//     Copyright notice: You are free to modify, extend and distribute
//        this code granted that the author of the original code is
//        mentioned as the original author of the code.

#include "linearAlgebra.h"
#include "mex.h"
#include "string.h"
#include "math.h"


void bubblesort(double *eta, double *xl, double *Pl, int size);

void fov_bubblesort(int *fov, double *eta, double *xl, double *Pl, int size);

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /*
     * nrhs     Number of input (right-side) arguments, or the size of the prhs array.
     * prhs     Array of input arguments.
     */
    
    /* variable declarations of input variables  */
    mxArray *xlPtr, *PlPtr, *etaPtr, *etaThresholdPtr;    /* field pointer */
    double *xn, *eta, *xl, *Pl;             /* input vectors and matrices */
    
    plhs[0] = mxDuplicateArray(prhs[0]);     // P_D, SP r, gating_size, lambda_c
    
    xlPtr = mxGetField(plhs[0], 0, "xl");
    PlPtr = mxGetField(plhs[0], 0, "Pl");
    etaPtr = mxGetField(plhs[0], 0, "eta");
    
    xn = mxGetPr(mxGetField(plhs[0], 0, "xn"));
    xl = mxGetPr(xlPtr);
    Pl = mxGetPr(PlPtr);
    eta = mxGetPr(etaPtr);
    
    double w_min = *mxGetPr(mxGetField(prhs[1], 0, "w_min"));
    double merging_threshold = *mxGetPr(mxGetField(prhs[1], 0, "merging_threshold"));
    int M_cap = (int)*mxGetPr(mxGetField(prhs[1], 0, "M_cap"));
    double FOV_RANGE2 = pow(*mxGetPr(mxGetField(prhs[1], 0, "fov_range")) +
            *mxGetPr(mxGetField(prhs[1], 0, "range_buffer")),2);
    double FOV_ANGLE = *mxGetPr(mxGetField(prhs[1], 0, "fov_angle")) +
            *mxGetPr(mxGetField(prhs[1], 0, "angle_buffer"));
    double etaT = log(mxGetScalar(mxGetField(prhs[1], 0, "P_B")));
    
    int xl_dim = (int)mxGetM(xlPtr);
    int n_k = (int)mxGetN(xlPtr);
    int xl_dim2 = xl_dim*xl_dim;
    
    // sort arrays based on the weight
    bubblesort(eta, xl, Pl, n_k);
    
    /* 
     *
     *  Hypothesis reduction, pruning
     *
     */
    int endofarray = n_k;
    for(int i = n_k; i > 0; --i)
    {
        if (eta[i] < w_min)
        {
            endofarray = i;
        }
    }
    
    double Sigma[4], detP[1], nu[2], x[2], P[4], Snu[2];
    double log_sum_w, d;

    char *merged_idx = mxCalloc(endofarray*sizeof(char),sizeof(char));
    char *cluster_idx = mxCalloc(endofarray*sizeof(char),sizeof(char));
    double *log_w = mxCalloc(endofarray*sizeof(double),sizeof(double));
    int *idx = mxCalloc(endofarray*sizeof(int),sizeof(int));
    int source_i, n, l, el = 0;
    
    /* 
     *
     *  Hypothesis reduction, merging
     *
     */
    for(int i = 0; i < endofarray; i++)
    {
        if(merged_idx[i] == 1)
            continue;
        
        cluster_idx[i] = 1;
        merged_idx[i] = 1;
        n = 1;
        
        // determine landmarks to be merged
        matrixInverse(Sigma, detP, Pl+(i*xl_dim2), xl_dim);
        for(int j = i+1; j < endofarray; j++)
        {
            nu[0] = xl[i*xl_dim+0] - xl[j*xl_dim+0];
            nu[1] = xl[i*xl_dim+1] - xl[j*xl_dim+1];
            
            matrixMultiply(xl_dim, xl_dim, xl_dim, 1, 0, 0, Sigma, nu, Snu);
            d = nu[0]*Snu[0] + nu[1]*Snu[1];
            
            if (d < merging_threshold)
            {
                cluster_idx[j] = 1;
                merged_idx[j] = 1;
                n++;
            }
        }
        
        // normalize log-weights
        log_sum_w = 0.0, l = 0;

        for(int j = i; j < endofarray; j++)
        {
            if(cluster_idx[j] == 1)
            {
                log_w[l] = eta[j];
                idx[l] = j;
                log_sum_w += exp(log_w[l]-log_w[0]);
                l++;
                cluster_idx[j] = 0;
            }
        }
        log_sum_w = log_w[0] + log(log_sum_w);
        for(int j = 0; j < n; j++)
            log_w[j] -= log_sum_w;
        
        
        // Moment matching
        x[0] = 0.0; x[1] = 0.0;
        P[0] = 0.0; P[1] = 0.0; P[2] = 0.0; P[3] = 0.0;
        
        // mean
        for(int j = 0; j < n; j++)
        {
            log_w[j] = exp(log_w[j]);
            for(int k = 0; k < xl_dim; k++)
                x[k] += log_w[j]*xl[idx[j]*xl_dim+k];
        }
        
        // covariance
        for(int j = 0; j < n; j++)
        {
            for(int jj = 0; jj < xl_dim; jj++)
                nu[jj] = x[jj] - xl[idx[j]*xl_dim+jj];
            
            for(int jj = 0; jj < xl_dim; jj++){
                for(int jjj = 0; jjj < xl_dim; jjj++)
                    P[jj*xl_dim+jjj] += log_w[j]*Pl[idx[j]*xl_dim2+jj*xl_dim+jjj] + log_w[j]*nu[jj]*nu[jjj];
            }
        }
        
        // covariance is very large, don't merge
        if (P[0] + P[3] > merging_threshold)
            log_sum_w = -41.4465;


        // store variables
        eta[el] = log_sum_w;
        memcpy(xl+xl_dim*el, x,xl_dim*sizeof(double));
        memcpy(Pl+xl_dim2*el,P,xl_dim2*sizeof(double));
        el++;
    }
     
    mxFree(merged_idx);
    mxFree(cluster_idx);
    mxFree(log_w);
    mxFree(idx);
    
    
    
    /*
     *
     *  Map management
     *
     */
    
    /* variable declarations of input variables  */
    mxArray *xlpPtr, *PlpPtr, *etapPtr, *etaTpPtr;    /* field pointer */
    double *etap, *etaTp, *xlp, *Plp;             /* input vectors and matrices */
    
    xlpPtr = mxGetField(plhs[0], 0, "xl_p");
    PlpPtr = mxGetField(plhs[0], 0, "Pl_p");
    etapPtr = mxGetField(plhs[0], 0, "eta_p");
    
    xlp = mxGetPr(xlpPtr);
    Plp = mxGetPr(PlpPtr);
    etap = mxGetPr(etapPtr);
    
    int np_k = (int)mxGetN(xlpPtr);
    
    int *in_fov = mxCalloc(el*sizeof(int),sizeof(int));
    int *out_fov = mxCalloc(np_k*sizeof(int),sizeof(int));
    int move_outside_fov = 0;
    int move_inside_fov = 0;
    double dx, dy;
    // determine active components that are inside/outside FOV
    for (int i = 0; i < el; i++)
    {
        dx = xl[i*2] - xn[0];
        dy = xl[i*2+1] - xn[1];
        if (pow(dx,2) + pow(dy,2) > FOV_RANGE2)
            in_fov[i] = 0;
        else if (fabs(pi_to_pi(atan2(dy,dx) - xn[2])) > FOV_ANGLE)
            in_fov[i] = 0;
        else
            in_fov[i] = 1;
        
        if (in_fov[i] == 0 && eta[i] <= etaT)
            in_fov[i] = -1;

        if (in_fov[i] == 0)
            move_outside_fov++;
    }
    
    // determine passive components that are inside/outside FOV
    for (int i = 0; i < np_k; i++)
    {
        dx = xlp[i*2] - xn[0];
        dy = xlp[i*2+1] - xn[1];
        if (pow(dx,2) + pow(dy,2) > FOV_RANGE2)
            out_fov[i] = 1;
        else if (fabs(pi_to_pi(atan2(dy,dx) - xn[2])) > FOV_ANGLE)
            out_fov[i] = 1;
        else
            out_fov[i] = 0;
        if (out_fov[i] == 0)
            move_inside_fov++;
    }
    
    
    fov_bubblesort(in_fov, eta, xl, Pl, el);
    int xl_allocated = 0;
    for (int i = 0; i < el; i++)
    {
        if(in_fov[i] == 1)
            xl_allocated++;
    }
    
    fov_bubblesort(out_fov, etap, xlp, Plp, np_k);
    int xlp_allocated = 0;
    for (int i = 0; i < np_k; i++)
    {
        if(out_fov[i] == 1)
            xlp_allocated++;
    }

    if(n_k-xl_allocated < move_inside_fov)
    {
        int xl_new_size = xl_allocated + move_inside_fov;
        
        // reallocate size of arrays
        mxSetData(xlPtr, (mxArray **)mxRealloc(mxGetData(xlPtr), xl_dim * xl_new_size * sizeof(mxArray *)));
        mxSetData(PlPtr, (mxArray **)mxRealloc(mxGetData(PlPtr), xl_dim2 * xl_new_size * sizeof(mxArray *)));
        mxSetData(etaPtr, (mxArray **)mxRealloc(mxGetData(etaPtr), xl_new_size * sizeof(mxArray *)));
        
        // get pointer to reallocated arrays
        xl = mxGetPr(xlPtr);
        Pl = mxGetPr(PlPtr);
        eta = mxGetPr(etaPtr);
    }
    
    if(np_k-xlp_allocated < move_outside_fov)
    {
        int xlp_new_size = xlp_allocated + move_outside_fov;
        
        if(xlp_allocated==0)
        {
            const mwSize ndim_init = 3;
            const mwSize dims_init[3] = {xl_dim,xl_dim,xlp_new_size};
            
            xlpPtr = mxCreateDoubleMatrix(xl_dim, xlp_new_size, mxREAL);
            PlpPtr = mxCreateNumericArray(ndim_init, dims_init, mxDOUBLE_CLASS, mxREAL);
            etapPtr = mxCreateDoubleMatrix(1, xlp_new_size, mxREAL);
            
            mxSetField(plhs[0], 0, "xl_p", xlpPtr);
            mxSetField(plhs[0], 0, "Pl_p", PlpPtr);
            mxSetField(plhs[0], 0, "eta_p", etapPtr);
        }
        else
        {
            // reallocate size of arrays
            mxSetData(xlpPtr, (mxArray **)mxRealloc(mxGetData(xlpPtr), xl_dim * xlp_new_size * sizeof(mxArray *)));
            mxSetData(PlpPtr, (mxArray **)mxRealloc(mxGetData(PlpPtr), xl_dim2 * xlp_new_size * sizeof(mxArray *)));
            mxSetData(etapPtr, (mxArray **)mxRealloc(mxGetData(etapPtr), xlp_new_size * sizeof(mxArray *)));
        }
        // get pointer to reallocated arrays
        xlp = mxGetPr(xlpPtr);
        Plp = mxGetPr(PlpPtr);
        etap = mxGetPr(etapPtr);
    }
    
    double *xl_tmp, *Pl_tmp, *eta_tmp;
    // make temporary copy of arrays
    
    xl += xl_dim * xl_allocated;
    Pl += xl_dim2 * xl_allocated;
    eta += xl_allocated;
    
    xlp += xl_dim * xlp_allocated;
    Plp += xl_dim2 * xlp_allocated;
    etap += xlp_allocated;
    
    if (move_outside_fov > 0)
    {
        xl_tmp = mxCalloc(move_outside_fov*xl_dim*sizeof(double),sizeof(double));
        Pl_tmp = mxCalloc(move_outside_fov*xl_dim2*sizeof(double),sizeof(double));
        eta_tmp = mxCalloc(move_outside_fov*sizeof(double),sizeof(double));
        
        memcpy(xl_tmp, xl,(xl_dim*move_outside_fov)*sizeof(double));
        memcpy(Pl_tmp, Pl,(xl_dim2*move_outside_fov)*sizeof(double));
        memcpy(eta_tmp, eta,move_outside_fov*sizeof(double));
    }
    
    if (move_inside_fov > 0)
    {
        memcpy(xl, xlp,(xl_dim*move_inside_fov)*sizeof(double));
        memcpy(Pl, Plp,(xl_dim2*move_inside_fov)*sizeof(double));
        memcpy(eta, etap,move_inside_fov*sizeof(double));
        
        xl_allocated += move_inside_fov;
    }
    
    if (move_outside_fov > 0)
    {
        memcpy(xlp, xl_tmp,(xl_dim*move_outside_fov)*sizeof(double));
        memcpy(Plp, Pl_tmp,(xl_dim2*move_outside_fov)*sizeof(double));
        memcpy(etap, eta_tmp,move_outside_fov*sizeof(double));
        
        xlp_allocated += move_outside_fov;
        
        mxFree(xl_tmp);
        mxFree(Pl_tmp);
        mxFree(eta_tmp);
    }
    
    const mwSize ndim = 3;
    const mwSize dims_xl[3] = {xl_dim,xl_dim,xl_allocated};
    const mwSize dims_xlp[3] = {xl_dim,xl_dim,xlp_allocated};
    
    // set the dimensions correctly
    mxSetN(xlPtr, xl_allocated);
    mxSetDimensions(PlPtr, dims_xl, ndim);
    mxSetN(etaPtr, xl_allocated);
    
    // set the dimensions correctly
    mxSetN(xlpPtr, xlp_allocated);
    mxSetDimensions(PlpPtr, dims_xlp, ndim);
    mxSetN(etapPtr, xlp_allocated);
}

void bubblesort(double *eta, double *xl, double *Pl, int size)
{
    double temp;
    int xl_dim = 2;
    int xl_dim2 = 4;
    
    // loop to access each array element
    for (int step = 0; step < size - 1; ++step) {
        
        // loop to compare array elements
        for (int i = 0; i < size - step - 1; ++i) {
            
            // compare two adjacent elements
            // change > to < to sort in descending order
            if (eta[i] < eta[i + 1]) {
                
                // swapping occurs if elements
                // are not in the intended order
                temp = eta[i];
                eta[i] = eta[i + 1];
                eta[i + 1] = temp;
                
                for (int j = 0; j < xl_dim; j++)
                {
                    temp = xl[i*xl_dim+j];
                    xl[i*xl_dim+j] = xl[(i + 1)*xl_dim+j];
                    xl[(i + 1)*xl_dim+j] = temp;
                }
                
                for (int j = 0; j < xl_dim2; j++)
                {
                    temp = Pl[i*xl_dim2+j];
                    Pl[i*xl_dim2+j] = Pl[(i + 1)*xl_dim2+j];
                    Pl[(i + 1)*xl_dim2+j] = temp;
                }
            }
        }
    }
}

void fov_bubblesort(int *fov, double *eta, double *xl, double *Pl, int size)
{
    double temp;
    int xl_dim = 2;
    int xl_dim2 = 4;
    
    // loop to access each array element
    for (int step = 0; step < size - 1; ++step) {
        
        // loop to compare array elements
        for (int i = 0; i < size - step - 1; ++i) {
            
            // compare two adjacent elements
            // change > to < to sort in descending order
            if (fov[i] < fov[i + 1]) {
                
                // swapping occurs if elements
                // are not in the intended order
                temp = fov[i];
                fov[i] = fov[i + 1];
                fov[i + 1] = temp;
                
                temp = eta[i];
                eta[i] = eta[i + 1];
                eta[i + 1] = temp;
                
                for (int j = 0; j < xl_dim; j++)
                {
                    temp = xl[i*xl_dim+j];
                    xl[i*xl_dim+j] = xl[(i + 1)*xl_dim+j];
                    xl[(i + 1)*xl_dim+j] = temp;
                }
                
                for (int j = 0; j < xl_dim2; j++)
                {
                    temp = Pl[i*xl_dim2+j];
                    Pl[i*xl_dim2+j] = Pl[(i + 1)*xl_dim2+j];
                    Pl[(i + 1)*xl_dim2+j] = temp;
                }
            }
        }
    }
}
