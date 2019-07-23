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
#include "includes\utilities.cpp"   // transformation and utility functions

// c. solve
#include "includes\interpolate.cpp"     // setup and interpolant gatewys
#include "includes\EGM.cpp"             // EGM


////////////////
// 4. gateway //
////////////////

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
     
        logs::solve(-1,"Solving model...\n"); // reset log file
        

    //////////////////
    // 5. setup par //
    //////////////////
    
    par_struct* par = new par_struct;
    par::setup(par,plhs,prhs,1);

        // index matrices
        auto indices = new int*[par->T];
        for(int t = 0; t < par->T; t++){

            indices[t] = new int[par->Nd[t]*par->NN*par->NP*3];
            for(int d   = 0; d   < par->Nd[t]; d++){
            for(int i_N = 0; i_N < par->NN;    i_N++){
            for(int i_P = 0; i_P < par->NP;    i_P++){

                int index = d*par->NN*par->NP*3 + i_N*par->NP*3 + i_P*3;
                indices[t][index + 0] = d;
                indices[t][index + 1] = i_N;
                indices[t][index + 2] = i_P;

            } } }

        }

    /////////////////
    // 6. parallel //
    /////////////////

    #pragma omp parallel num_threads(THREADS)
    {

        // a. setup workers
        auto sol = new sol_struct;
        sol::setup(par,sol);

        // b. time loop
        for(int t = par->T-1; t >= par->tmin; t--){

                #pragma omp master
                logs::solve(0,"t = %d\n", t);

            sol->t = t;

            // i. create interpolants
            if(t == par->T-2){
                interpolate::setup(t+1,par,sol,par->NA);
            } else if(t < par->T-2){
                interpolate::update(t+1,par,sol);                
            }

            // ii. precomputations
            EGM::all(t,par,sol,indices[t]);
  
        }

        // c. clean up workers
        interpolate::destroy(par,sol);        
        sol::destroy(sol);
        delete sol;

    } // parallel


    /////////////////
    // 7. clean up //
    /////////////////
    
    for(int t = 0; t < par->T; t++){
        delete[] indices[t]; 
    }
    delete[] indices;

    par::destroy(par,1);
    delete par;


        logs::solve(0,"Done.\n");

        // clean assertions file
        FILE* log_file = fopen("log_assert.txt","w");
        fprintf(log_file,"\n");
        fclose(log_file);

} // mex gateway