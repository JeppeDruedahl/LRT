namespace interpolate {

//////////////
// 1. setup //
//////////////

void setup(int t, par_struct *par, sol_struct *sol, int Nxi){

	int level = 3;

	for(int d   = 0; d   < par->Nd_max; d++){   // discrete states
	for(int i_P = 0; i_P < par->NP;    	i_P++){
	for(int i_N = 0; i_N < par->NN;    	i_N++){

		int i_cell = index::d2(t,d,par->Nd_max);
		int i_grid = index::d3(i_P,i_N,0,par->NN,par->NM);

		// a. grids				
		int dimx = 1;
		int Nx[1]; Nx[0] = par->NM;
		double *x[1]; x[0] = &par->M_ast[i_cell][i_grid];
		
		// b. values
		int dimy = 2;

			double *y[2];
			y[0] = &par->C_ast[i_cell][i_grid];
			y[1] = &par->V_ast[i_cell][i_grid];

		// c. create
		int i_interp = index::d3(d,i_P,i_N,par->NP,par->NN);			
		sol->interp[i_interp] = new linear_interp::set_struct;
		auto interp = sol->interp[i_interp];
        linear_interp::create(interp, dimx, dimy, Nx, x, y, Nxi);
			
			if(Nxi > 1) { // only print when solving
				logs::solve(level+1,"create interpolation (d=%d,i_P=%d,i_N=%d)\n",d,i_P,i_N);
			}
			
	} } }

}

void update(int t, par_struct *par, sol_struct *sol){

	int level = 3;

	for(int d   = 0; d   < par->Nd_max; d++){   // discrete states
	for(int i_P = 0; i_P < par->NP;    	i_P++){
	for(int i_N = 0; i_N < par->NN;    	i_N++){

		int i_cell = index::d2(t,d,par->Nd_max);
		int i_grid = index::d3(i_P,i_N,0,par->NN,par->NM);

		int i_interp = index::d3(d,i_P,i_N,par->NP,par->NN);			
		auto interp = sol->interp[i_interp];

		interp->Nx[0] = par->NM;
        interp->x[0]  = &par->M_ast[i_cell][i_grid];
   
	    interp->y[0] = &par->C_ast[i_cell][i_grid];
	    interp->y[1] = &par->V_ast[i_cell][i_grid];

	} } }

}


void destroy(par_struct *par, sol_struct *sol){

	for(int d   = 0; d   < par->Nd_max; d++){   // discrete states
	for(int i_P = 0; i_P < par->NP;    	i_P++){
	for(int i_N = 0; i_N < par->NN;    	i_N++){

		int i_interp = index::d3(d,i_P,i_N,par->NP,par->NN);
		linear_interp::destroy(sol->interp[i_interp]);
		
	} } }

}

/////////////////////
// 2. interpolants //
/////////////////////

void evaluate(sol_struct *sol, int t, int d, 
			  double P, double N, double *M,
 			  double *C, double *V, int Nxi){

	int level = 6;

		logs::solve(level,"interpolate::evaluate (d = %d, P = %g, N = %g, M[0] = %g)\n",d,P,N,M[0]);

	// a. unpack
	auto par = sol->par;
	
	auto P_vec = &par->grid_P[index::d2(t,0,par->NP)];
	auto N_vec = &par->grid_N[index::d2(t,0,par->NN)];	
	
	auto C_tmp = sol->C_tmp;
	auto V_tmp = sol->V_tmp;

	// b. xi
	double xi[1];
	xi[0] = M[0];
	auto xi_vec = M;

	// c. find Ns
	int i_N_left = linear_interp::binary_search(0, par->NN, N_vec, N);
	double N_reldiff = (N - N_vec[i_N_left])/(N_vec[i_N_left+1]-N_vec[i_N_left]);

		logs::solve(level,"N_reldiff = %g [%g %g]\n",N_reldiff,N_vec[i_N_left],N_vec[i_N_left+1]);

	// d. LRT interpolation
	if(par->LRT == 1){

		// i. interpolation in M
		for(int add = 0; add < 2; add++){

			int i_N = i_N_left + add;

			int i_interp = index::d3(d,0,i_N,par->NP,par->NN);
		    auto interp = sol->interp[i_interp];

			linear_interp::evaluate(interp,&sol->yi[add*2*Nxi],xi,xi_vec);

		}

	    // ii. interpolation in N
	    double *C_left 	= &sol->yi[0]; 
	    double *C_right = &sol->yi[2*Nxi];
	    double *V_left 	= &sol->yi[Nxi];
	    double *V_right = &sol->yi[3*Nxi];
	    	 
	    for(int i_A = 0; i_A < Nxi; i_A++){
	    	
	    	C[i_A] = C_left[i_A] + N_reldiff*(C_right[i_A]- C_left[i_A]);
	    	V[i_A] = V_left[i_A] + N_reldiff*(V_right[i_A]- V_left[i_A]);

			logs::solve(level+2,"M = %g, C = %g [%g %g], V = %g [%g %g]\n",M[i_A],C[i_A],C_left[i_A],C_right[i_A],
				V[i_A],V_left[i_A],V_right[i_A]);
				
	    }

		return;

	}

	// e. find Ps
	int i_P_left = linear_interp::binary_search(0, par->NP, P_vec, P);
	double P_reldiff = (P - P_vec[i_P_left])/(P_vec[i_P_left+1]-P_vec[i_P_left]);
	
		logs::solve(level,"P_reldiff = %g [%g %g]\n",P_reldiff,P_vec[i_P_left],P_vec[i_P_left+1]);

    // f. interpolation in (M,N)    
	for(int add_P = 0; add_P < 2; add_P++){

		int i_P = i_P_left + add_P;

		// i. interpolation in M
		for(int add_N = 0; add_N < 2; add_N++){

			int i_N = i_N_left + add_N;

			int i_interp = index::d3(d,i_P,i_N,par->NP,par->NN);
		    auto interp = sol->interp[i_interp];

			linear_interp::evaluate(interp,&sol->yi[add_N*2*Nxi],xi,xi_vec);		
			
		} 

		// ii. interpolation in N
	    double *C_left 	= &sol->yi[0]; 
	    double *C_right = &sol->yi[2*Nxi];
	    double *V_left 	= &sol->yi[Nxi];
	    double *V_right = &sol->yi[3*Nxi];
	    	 
	    for(int i_A = 0; i_A < Nxi; i_A++){
	    	C_tmp[add_P*Nxi + i_A] = C_left[i_A] + N_reldiff*(C_right[i_A]- C_left[i_A]);
	    	V_tmp[add_P*Nxi + i_A] = V_left[i_A] + N_reldiff*(V_right[i_A]- V_left[i_A]);
			logs::solve(level+2,"add_P = %d: M = %g, C = %g [%g %g], V = %g [%g %g]\n",add_P,M[i_A],C_tmp[add_P*Nxi + i_A],C_left[i_A],C_right[i_A],
				V_tmp[add_P*Nxi + i_A],V_left[i_A],V_right[i_A]);
		}

	
    }

    // g. interpolation in P
    double *C_left 	= &C_tmp[0]; 
    double *C_right = &C_tmp[Nxi];
    double *V_left 	= &V_tmp[0];
    double *V_right = &V_tmp[Nxi];

    for(int i_A = 0; i_A < Nxi; i_A++){
    	C[i_A] = C_left[i_A] + P_reldiff*(C_right[i_A] - C_left[i_A]);
    	V[i_A] = V_left[i_A] + P_reldiff*(V_right[i_A] - V_left[i_A]);
		logs::solve(level+1,"M = %g, C = %g [%g %g], V = %g [%g %g]\n",M[i_A],C[i_A],C_left[i_A],C_right[i_A],
			V[i_A],V_left[i_A],V_right[i_A]);
    }
	
}

} // namespace