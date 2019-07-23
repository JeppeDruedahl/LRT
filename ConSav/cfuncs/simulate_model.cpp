//////////////////////////
// 1. external includes //
//////////////////////////

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <string>
#include <omp.h>
#include "mex.h"
#include "matrix.h"


//////////////////////////
// 2. define statements //
//////////////////////////

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define BOUND(X,A,B) MIN(MAX(X,A),B)
#define THREADS MIN(MAXTHREADS, omp_get_max_threads()-1)


//////////////////////////
// 3. internal includes //
//////////////////////////

// a. generic
#include "includes\assert.cpp"             // assert() function
#include "includes\logs.cpp"               // log:: functions
#include "includes\linear_interp.cpp"      // linear_interp:: functions
#include "includes\index.cpp"              // index:: functions
#include "includes\mymex.cpp"              // functions to interact with mex

// b. basic
#include "includes\par_struct.cpp"  // define par_struct + setup/destroy functions
#include "includes\sol_struct.cpp"  // define sol_struct + setup/destroy functions

// c. solve
#include "includes\interpolate.cpp"     // setup and interpolant gatewys


////////////////
// 4. gateway //
////////////////

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    logs::simulate(-1,"Starting simulation...\n");


    //////////////////
    // 5. setup par //
    //////////////////
    
    par_struct* par = new par_struct;
    par::setup(par,plhs,prhs,2);

        logs::simulate(0,"Setup completed.\n");

        // easy access to outputs
        auto d = par->d;
        
        auto M = par->M;
        auto N = par->N;
        auto P = par->P;              
        auto Y = par->Y;    
        auto C = par->C;
        auto A = par->A; 
        auto V = par->V; 
        auto MPC = par->MPC; 

    /////////////////
    // 6. parallel //
    /////////////////
        
    #pragma omp parallel num_threads(THREADS)
    {

        auto sol = new sol_struct;
        sol::setup(par,sol);  


    //////////////////
    // 7. time loop //
    //////////////////

    for(int t = 0; t < par->simT; t++){ // forward through time

        #pragma omp master
        logs::simulate(0,"t = %d\n",t);
        
        // create interpolants evaluate for a single x
        if(par->sim_only_inc == 0){
            if(t == 0){
                interpolate::setup(t,par,sol,1);
            } else {
                interpolate::update(t,par,sol);
            }
        }

    /////////////////////////
    // 8. individuals loop //
    /////////////////////////

        logs::simulate(1,"individuals loop\n");

    #pragma omp for
    for(int i = 0; i < par->simN; i++){
	
            logs::simulate(2,"i = %d\n",i);

        // a. index
        int index      = t*par->simN + i;
        int index_plus = (t+1)*par->simN + i;
	
        // b. initial values
        if(t == 0){

            d[index] = par->d_ini[i];
            if(par->LRT == 1){
                P[index] = NAN;
                Y[index] = par->grid_Y[0*par->Nd_max+d[index]];
            } else {
                P[index] = par->P_ini[i]*par->psi_sim[i];
                Y[index] = P[index]*par->xi_sim[i];                        
            }
            M[index] = par->M_ini[i] + Y[index];
            N[index] = par->kappa*Y[index];

        }

            logs::simulate(2,"states: d = %d, M = %g, N = %g, P = %g\n",
                d[index],M[index],N[index],P[index]);

        // c. optimal choice
        if(par->sim_only_inc == 0){

            interpolate::evaluate(sol, t, d[index],
                                   P[index], N[index], &M[index],
                                   &C[index], &V[index], 1);     
            if(par->epstein_zin == 0){
                V[index] = -1.0/V[index];
            }

                // MPC
                if(t < par->simT){
                
                    double M_add[1], C_add[1], V_add[1];
                    double MPC_add = 0.01*Y[index];
                    M_add[0] = M[index] + MPC_add;
                    interpolate::evaluate(sol, t, d[index],
                                           P[index], N[index], M_add,
                                           C_add, V_add, 1);  
                    
                    MPC[index] = (C_add[0]-C[index])/MPC_add;
                
                } else {

                    MPC[index] = NAN;

                }

            A[index] = M[index] - C[index];

                logs::simulate(2,"choices: C = %g, A = %g, M = %g\n",C[index],A[index],M[index]);

        }

        // d. next period discrete state income
        if(t == par->simT-1){ continue; }
        if(par->LRT == 1){

            int i_d_states = index::d2(t,d[index],par->Nd_max);
            
            // i. next period groups
            int Nd_plus      = par->Nd_plus[i_d_states];
            int* i_d_plus    = par->i_d_plus[i_d_states];
            double* p_cumsum = par->i_d_plus_p_cumsum[i_d_states];
            double* grid_Y   = &par->grid_Y[(t+1)*par->Nd_max];
            
            // ii. transition
            double p = par->p_sim[index_plus];
            int j = 0;
            while(p > p_cumsum[j] && j+1 <= Nd_plus-1){
                j++;
            }

            // iii. group
            d[index_plus] = i_d_plus[j];

            // iv. income
            Y[index_plus] = grid_Y[d[index_plus]];

        } else {

            d[index_plus] = 0;
            P[index_plus] = par->Gamma[t]*P[index]*par->psi_sim[index_plus];
            Y[index_plus] = P[index_plus]*par->xi_sim[index_plus];
        
        }

            logs::simulate(2,"income: Y_plus = %g\n",Y[index_plus]);

        // e. next period continuous states
        M[index_plus] = par->R*A[index] + exp(par->inceq)*Y[index_plus];
        M[index_plus] = MAX(M[index_plus],1e-6);
    
        N[index_plus] = par->Rb*N[index] + par->kappa*Y[index_plus];
	
	
    } // N
    #pragma omp barrier

    } // T        
        
        if(par->sim_only_inc == 0){
            interpolate::destroy(par,sol);
        }

        sol::destroy(sol);
        delete sol;
	
    } // parallel


    //////////////////
    // 9. clean up //
    //////////////////

    par::destroy(par,2);
    delete par;

        logs::simulate(0,"Done.\n");

        // clean assertions file
        FILE* log_file = fopen("log_assert.txt","w");
        fprintf(log_file,"\n");
        fclose(log_file);

} // simulate