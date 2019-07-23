namespace EGM {

    // forward declarations
    void last_period(sol_struct *sol);
    void EGM(sol_struct *sol);
    void next_period_values(sol_struct *sol);
    void update_W_and_Q(sol_struct *sol);

////////////////
// 1. profile //
////////////////

void all(int t, par_struct *par, sol_struct *sol, int* indices){
	
    int level = 2;
        
        logs::solve(level,"EGM::all\n");

    // number of points at borrowing constraint
    int NBC = par->NM-par->NA;

    // parallel
    #pragma omp for schedule(dynamic)
    for(int i = 0; i < par->Nd[t]*par->NN*par->NP; i++){

        int d = indices[3*i + 0];
        int i_N = indices[3*i + 1];
        int i_P = indices[3*i + 2];

            logs::solve(level+1,"(d,i_P,i_N) = (%1d,%3d,%3d)\n",d,i_P,i_N);

    	// a. pack
    	sol->t = t;
    	sol->d = d;
    	sol->P = par->grid_P[index::d2(t,i_P,par->NP)];
        sol->N = par->grid_N[index::d2(t,i_N,par->NN)];

        // b. post-decision vector
        sol->A_vec = &par->grid_A[index::d2(t,0,par->NA)];

        // c. last-period
        if(t == par->T-1){

            int i_cell    = index::d2(t,d,par->Nd_max);        
            int i_states  = index::d3(i_P,i_N,0,par->NN,par->NM);
            
            sol->M_vec = &par->M_ast[i_cell][i_states];        
            sol->C_vec = &par->C_ast[i_cell][i_states];     
            sol->V_vec = &par->V_ast[i_cell][i_states];       

            last_period(sol);

            continue;

        } // last period

        // d. EGM
        if(par->LRT == 1){
        int i_d_states = index::d2(t,d,par->Nd_max);
        
            sol->Nd_plus    = par->Nd_plus[i_d_states];
            sol->i_d_plus   = par->i_d_plus[i_d_states];
            sol->i_d_plus_p = par->i_d_plus_p[i_d_states];
        
        }

        int i_cell    = index::d2(t,d,par->Nd_max);        
        int i_states  = index::d3(i_P,i_N,NBC,par->NN,par->NM);
        
            sol->M_vec = &par->M_ast[i_cell][i_states];        
            sol->C_vec = &par->C_ast[i_cell][i_states];     
            sol->V_vec = &par->V_ast[i_cell][i_states];        

        EGM(sol);

        // e. constraint
        i_states  = index::d3(i_P,i_N,0,par->NN,par->NM);
        
            sol->M_vec = &par->M_ast[i_cell][i_states];        
            sol->C_vec = &par->C_ast[i_cell][i_states];     
            sol->V_vec = &par->V_ast[i_cell][i_states];   

        for(int i_M = 0; i_M < NBC; i_M++){
            sol->M_vec[i_M] = par->Mmin + (double)i_M/double(NBC)*sol->M_vec[NBC];
            sol->C_vec[i_M] = sol->M_vec[i_M];
            sol->V_vec[i_M] = utilities::value_of_choice(sol->C_vec[i_M],sol->W_vec[0],par);
            if(par->epstein_zin == 0){
                sol->V_vec[i_M] = -1.0/sol->V_vec[i_M];
            }
        }

    }

}


void last_period(sol_struct *sol){

    int level = 3;

        logs::solve(level,"EGM::last_period\n");    

    //////////////
    // 1. setup //
    //////////////

    // a. unpack
    auto par   = sol->par;
    auto t     = sol->t;
    auto A_vec = sol->A_vec;
    auto C_vec = sol->C_vec;
    auto M_vec = sol->M_vec;
    auto V_vec = sol->V_vec; 

    double M_retire_max = par->grid_A[index::d2(par->T-1,par->NA-1,par->NA)];
    for(int i_M = 0; i_M < par->NM; i_M++){

        // a. cash-on-hand
        M_vec[i_M] = par->Mmin + (double)i_M/double(par->NM)*M_retire_max;
        
        // b. consumption
        if(par->kappa == 0.0){
            C_vec[i_M] = M_vec[i_M]/par->zeta;
        } else {
            C_vec[i_M] = (M_vec[i_M] + sol->N)/par->zeta;            
        }

        // c. value
        if(par->epstein_zin == 0){
            V_vec[i_M] = (1.0-par->beta)*pow(C_vec[i_M],1.0-par->rho)/(1.0-par->rho);
            V_vec[i_M] = -1.0/V_vec[i_M];            
        } else {
            V_vec[i_M] = pow((1.0-par->beta),1.0/(1.0-par->sigma))*C_vec[i_M];
        }
    
    }


}

////////////
// 2. EGM //
////////////

void EGM(sol_struct *sol){

    int level = 3;

        logs::solve(level,"EGM::EGM\n");    

    //////////////
    // 1. setup //
    //////////////

    // a. unpack
	auto par   = sol->par;
	auto t     = sol->t;
    auto A_vec = sol->A_vec;
    auto C_vec = sol->C_vec;
    auto M_vec = sol->M_vec;
    auto V_vec = sol->V_vec;    
    auto W_vec = sol->W_vec;
    auto Q_vec = sol->Q_vec;

    /////////////
    // 2. loop //
    /////////////
    
    // a. initialize
	for(int i_A = 0; i_A < par->NA; i_A++){
	   W_vec[i_A] = 0.0;
	   Q_vec[i_A] = 0.0;
	}
        
    // b. loop over shocks
    int Nshocks;
    if(par->LRT == 1){
        Nshocks      = sol->Nd_plus;
        sol->weights = sol->i_d_plus_p;
    } else {
    	Nshocks      = par->Nshocks[t+1];
        sol->weights = &par->weights[(t+1)*par->Nshocks_max];
        sol->psi     = &par->psi[(t+1)*par->Nshocks_max];       
        sol->xi      = &par->xi[(t+1)*par->Nshocks_max];
    }

	for(int i_shock = 0; i_shock < Nshocks; i_shock++){

        logs::solve(level+2,"i_shock = %d\n",i_shock);

        sol->i_shock = i_shock;

        // i. next-period states
    	next_period_values(sol);

        // ii. interpolaion
    	interpolate::evaluate(sol, t+1, sol->d_plus,
                              sol->P_plus, sol->N_plus, sol->M_plus,
                              sol->C_plus_vec, sol->V_plus_vec, par->NA);

        // iii. update
        update_W_and_Q(sol);

    }

    /////////////////////////
    // 3. finalize W and Q //
    /////////////////////////

    if(par->epstein_zin == 1){
        for(int i_A = 0; i_A < par->NA; i_A++){
            W_vec[i_A] = pow(W_vec[i_A],1.0/(1.0-par->rho));            
            Q_vec[i_A] = Q_vec[i_A]*pow(W_vec[i_A],par->rho-par->sigma);
        }
    }


    ///////////////////
    // 4. C, M and V //
    ///////////////////

    for(int i_A = 0; i_A < par->NA; i_A++){

        // a. invert Euler-equation and endogenous M
        C_vec[i_A] = utilities::inv_marg_u(Q_vec[i_A],par);            
            
        // b. endogenous cash-on-hand
        M_vec[i_A] = A_vec[i_A] + C_vec[i_A];

        // c. value function
        V_vec[i_A] = utilities::value_of_choice(C_vec[i_A],W_vec[i_A],par);
		if(par->epstein_zin == 0){
            V_vec[i_A] = -1.0/V_vec[i_A];
        }	

        logs::solve(level+2,"A = %g, Q = %g, W = %d, C = %g, V = %\n",
            A_vec[i_A],Q_vec[i_A],W_vec[i_A],C_vec[i_A],V_vec[i_A]);

    }

}

void next_period_values(sol_struct *sol){

    int level = 5;

        logs::solve(level,"EGM::next_period_values\n");

    // a. unpack
    auto par = sol->par;
    auto t   = sol->t;

    // b. income
    if(par->LRT == 1){
    
        sol->P_plus  = NAN;
        sol->d_plus  = sol->i_d_plus[sol->i_shock];
        sol->Y_plus  = par->grid_Y[(t+1)*par->Nd_max + sol->d_plus];

    } else {

		sol->d_plus = 0;
	
        double psi = sol->psi[sol->i_shock];
        double xi  = sol->xi[sol->i_shock];

        double P_min    = par->grid_P[(t+1)*par->NP];
        double P_max    = par->grid_P[(t+1)*par->NP + par->NP-1];

        double P_plus   = par->Gamma[t]*sol->P*psi;
        sol->P_plus     = BOUND(P_plus,P_min,P_max);
        sol->Y_plus     = sol->P_plus*xi;

    }

    // c. illiquid assets
    sol->N_plus = par->Rb*sol->N + par->kappa*sol->Y_plus;

    // d. cash-on-hand
    for(int i_A = 0; i_A < par->NA; i_A++){

        sol->M_plus[i_A] = par->R*sol->A_vec[i_A] + exp(par->inceq)*sol->Y_plus;
        sol->M_plus[i_A] = MAX(sol->M_plus[i_A],1e-6);
		
		  logs::solve(level+1,"A = %g, Y_plus = %g, M_plus = %g\n",sol->A_vec[i_A],sol->Y_plus,sol->M_plus[i_A]);
		
    }
    
}

void update_W_and_Q(sol_struct *sol){

    int level = 5;

        logs::solve(level,"EGM::update_W_and_Q\n");

    auto par = sol->par;

    for(int i_A = 0; i_A < par->NA; i_A++){
		
		double W_shocks, Q_shocks;
        double weight = sol->weights[sol->i_shock];
 
        // 1. update post-decision value function
        double V_plus = sol->V_plus_vec[i_A];
        if(par->epstein_zin == 0){ // CRRA
            W_shocks = -1.0/V_plus;
        } else { // Epstein-Zin
            W_shocks = pow(V_plus,1.0-par->rho); 
        }
            
        // 2. update post-decision marginal utility of cash
        double C_plus = sol->C_plus_vec[i_A];
        if(par->epstein_zin == 0){ // CRRA
            Q_shocks = par->beta*par->R*utilities::marg_u(C_plus,par);            
        } else { // Epstein-Zin
            Q_shocks = par->beta*par->R*utilities::marg_u(C_plus,par)*pow(V_plus,par->sigma-par->rho); 
        } 
		
		sol->W_vec[i_A] += weight*W_shocks;
		sol->Q_vec[i_A] += weight*Q_shocks;
		
		logs::solve(level+1,"weight = %g, W_shocks = %g, Q_shocks = %g\n",weight,W_shocks,Q_shocks);
		
    }

}


} // namespace