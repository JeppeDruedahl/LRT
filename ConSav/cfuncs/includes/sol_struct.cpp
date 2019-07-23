//////////////////
// 1. variables //
//////////////////

typedef struct
{
 
    // a. par
    par_struct *par;

    // b. interpolants
    linear_interp::set_struct **interp;
    double *yi, *C_tmp, *V_tmp;

    // c. info
    int t, d, d_plus;
    double P, N, Y_plus, P_plus, N_plus, *M_plus;
    double *psi, *xi, *weights;   
    
    // d. EGM
    int i_shock, Nd_plus, *i_d_plus;
    double *i_d_plus_p; 
    double *A_vec, *C_plus_vec, *V_plus_vec;
    double *C_vec, *M_vec, *V_vec;
    double *W_vec, *Q_vec;

} sol_struct;


namespace sol {

//////////////
// 2. setup //
//////////////

void setup(par_struct *par, sol_struct *sol){

    // a. par
    sol->par = par;
          
     // b. interpolants
    sol->interp = new linear_interp::set_struct*[par->Nd_max*par->NP*par->NN];

        sol->yi = new double[2*2*par->NA];
        sol->C_tmp = new double[2*par->NA];
        sol->V_tmp = new double[2*par->NA];        

    // c. EGM

        sol->M_plus = new double[par->NA];
        sol->C_plus_vec = new double[par->NA];
        sol->V_plus_vec = new double[par->NA];
        sol->W_vec = new double[par->NA];
		sol->Q_vec = new double[par->NA];

}


////////////////
// 3. destroy //
////////////////

void destroy(sol_struct *sol){

    // b. interpolants
    delete[] sol->interp;
    delete[] sol->yi;  
    delete[] sol->C_tmp;
    delete[] sol->V_tmp;

    // c. EGM
    delete[] sol->M_plus;
    delete[] sol->C_plus_vec;
    delete[] sol->V_plus_vec;    
    delete[] sol->W_vec;
    delete[] sol->Q_vec;

}

} // namespace