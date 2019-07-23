//////////////////
// 1. variables //
//////////////////

typedef struct{

    // a. grids
    int NA, NM, NN, NP, Nd_max, *Nd;  
    double *grid_A, *grid_N, *grid_P, *grid_Y;
    double Mmin;

    // b. demographics
    int T;

    // c. preferences
    int epstein_zin;
    double beta, rho, sigma, zeta;

    // d. new income process  
    int *Nd_plus;
    int **i_d_plus;
    double **i_d_plus_p, **i_d_plus_p_cumsum;

    // e. old income process      
    double *Gamma;        
    int Nshocks_max, *Nshocks;
    double *psi, *xi, *weights;    

    // f. assets
    double R, Rb;
    double kappa;

    // g. technical
    int LRT, tmin, sim_only_inc;
    double inceq;

    // h. simulate
    int simT, simN;
    int *d_ini;
    double *P_ini, *M_ini, *N_ini, *Y_ini;
    double *p_sim, *psi_sim, *xi_sim;

    // output
        
        // solve
        double **C_ast, **M_ast, **V_ast;

        // simulate
        int *d;
        double *M, *N, *Y, *P, *C, *A, *V, *MPC;

} par_struct;


namespace par {

//////////////
// 2. setup //
//////////////

void setup(par_struct *par, mxArray *plhs[], const mxArray *prhs[], int type){

    ///////////////
    // 1. inputs //
    ///////////////

    // a. grids   
    par->NA     = (int) mxGetScalar(mxGetField(prhs[0],0,"NA"));
    par->NM     = (int) mxGetScalar(mxGetField(prhs[0],0,"NM"));
    par->NN     = (int) mxGetScalar(mxGetField(prhs[0],0,"NN"));  
    par->NP     = (int) mxGetScalar(mxGetField(prhs[0],0,"NP"));             
    par->Nd_max = (int) mxGetScalar(mxGetField(prhs[0],0,"Nd_max"));
    par->Nd     = (int*) mxGetData(mxGetField(prhs[0],0,"Nd"));     
    par->grid_A = (double*) mxGetPr(mxGetField(prhs[0],0,"grid_A"));
    par->grid_N = (double*) mxGetPr(mxGetField(prhs[0],0,"grid_N"));    
    par->grid_P = (double*) mxGetPr(mxGetField(prhs[0],0,"grid_P"));
    par->grid_Y = (double*) mxGetPr(mxGetField(prhs[0],0,"grid_Y"));
    par->Mmin   = (double) mxGetScalar(mxGetField(prhs[0],0,"Mmin"));

    // b. demographics
    par->T   = (int) mxGetScalar(mxGetField(prhs[0],0,"T"));

    // c. preferences
    par->epstein_zin  = (int) mxGetScalar(mxGetField(prhs[0],0,"epstein_zin"));

    par->beta  = (double) mxGetScalar(mxGetField(prhs[0],0,"beta"));
    par->rho   = (double) mxGetScalar(mxGetField(prhs[0],0,"rho"));
    par->sigma = (double) mxGetScalar(mxGetField(prhs[0],0,"sigma"));
    par->zeta = (double) mxGetScalar(mxGetField(prhs[0],0,"zeta"));

    // d. new income process
    par->Nd_plus = (int*) mxGetData(mxGetField(prhs[0],0,"Nd_plus"));       

    // e. old income process
    par->Gamma       = (double*) mxGetPr(mxGetField(prhs[0],0,"Gamma"));
    par->Nshocks_max = (int) mxGetScalar(mxGetField(prhs[0],0,"Nshocks_max"));
    par->Nshocks     = (int*) mxGetData(mxGetField(prhs[0],0,"Nshocks"));
    par->psi         = (double*) mxGetPr(mxGetField(prhs[0],0,"psi"));
    par->xi          = (double*) mxGetPr(mxGetField(prhs[0],0,"xi"));      
    par->weights     = (double*) mxGetPr(mxGetField(prhs[0],0,"weights"));

    // f. assets
    par->R = (double) mxGetScalar(mxGetField(prhs[0],0,"R"));
    par->Rb = (double) mxGetScalar(mxGetField(prhs[0],0,"Rb"));    
    par->kappa = (double) mxGetScalar(mxGetField(prhs[0],0,"kappa"));

    // g. technical
    par->LRT = (int) mxGetScalar(mxGetField(prhs[0],0,"LRT"));
    par->tmin = (int) mxGetScalar(mxGetField(prhs[0],0,"tmin"));
    par->inceq = (double) mxGetScalar(mxGetField(prhs[0],0,"inceq"));
    par->sim_only_inc = (int) mxGetScalar(mxGetField(prhs[0],0,"sim_only_inc")); 

    // h. simulate
    if(type == 2){
    
        par->simN = (int) mxGetScalar(mxGetField(prhs[0],0,"simN"));
        par->simT = (int) mxGetScalar(mxGetField(prhs[0],0,"simT"));        
        par->d_ini = (int*) mxGetData(mxGetField(prhs[0],0,"d_ini"));
        par->P_ini = (double*) mxGetPr(mxGetField(prhs[0],0,"P_ini"));
        par->Y_ini = (double*) mxGetPr(mxGetField(prhs[0],0,"Y_ini"));        
        par->M_ini = (double*) mxGetPr(mxGetField(prhs[0],0,"M_ini"));
        par->N_ini = (double*) mxGetPr(mxGetField(prhs[0],0,"N_ini"));                
        par->p_sim = (double*) mxGetPr(mxGetField(prhs[0],0,"p_sim"));
        par->psi_sim = (double*) mxGetPr(mxGetField(prhs[0],0,"psi_sim"));
        par->xi_sim = (double*) mxGetPr(mxGetField(prhs[0],0,"xi_sim"));
    
    }

        if(type == 1){
            logs::solve(0," scalars and pointers loaded\n");        
        } else {
            logs::simulate(0," scalars and pointers loaded\n");                    
        }

    ////////////////////
    // 2. cell inputs //
    ////////////////////
    
    if(par->LRT == 1){

        par->i_d_plus          = new int*[par->Nd_max*par->T];
        par->i_d_plus_p        = new double*[par->Nd_max*par->T];
        par->i_d_plus_p_cumsum = new double*[par->Nd_max*par->T];

        auto i_d_plus_cell           = mxGetField(prhs[0],0,"i_d_plus");
        auto i_d_plus_p_cell         = mxGetField(prhs[0],0,"i_d_plus_p");
        auto i_d_plus_p_cumsum_cell  = mxGetField(prhs[0],0,"i_d_plus_p_cumsum");

        for(int t = 0; t < par->T-1;   t++){
        for(int j = 0; j < par->Nd[t]; j++){        

            int index = t*par->Nd_max + j;
        
            par->i_d_plus[index]           = (int*) mxGetData(mxGetCell(i_d_plus_cell,index));
            par->i_d_plus_p[index]         = (double*) mxGetPr(mxGetCell(i_d_plus_p_cell,index));
            par->i_d_plus_p_cumsum[index]  = (double*) mxGetPr(mxGetCell(i_d_plus_p_cumsum_cell,index));

        } }

        if(type == 1){        
            logs::solve(0," cells loaded\n");        
        } else {
            logs::simulate(0," cells loaded\n");                    
        }

    }

    ////////////////////////
    // 3. outputs - solve //
    ////////////////////////

    if(type == 1){

        // a. struct
        const char *field_names[] = {"C_ast", "M_ast", "V_ast"};
        
        int num_fields = sizeof(field_names)/sizeof(*field_names);
        plhs[0] = mxCreateStructMatrix(1, 1, num_fields, field_names);
        auto sol_struct = plhs[0];

        // b. cell dimensions
        int ndim_cell  = 2;
        auto dims_cell = new size_t[2];

        // c. array dimensions
        auto ndim = new size_t[1];
        auto dims = new size_t*[1];
                        
        // d. solution

            // cell
            dims_cell[0] = par->Nd_max;
            dims_cell[1] = par->T;

            // array      
            ndim[0] = 3;
            dims[0] = new size_t[3];           
            dims[0][0] = par->NM;
            dims[0][1] = par->NN;
            dims[0][2] = par->NP;

        par->C_ast = mymex::set_field_cell(sol_struct,"C_ast",ndim_cell,dims_cell,ndim,dims); 
        par->M_ast = mymex::set_field_cell(sol_struct,"M_ast",ndim_cell,dims_cell,ndim,dims);
        par->V_ast = mymex::set_field_cell(sol_struct,"V_ast",ndim_cell,dims_cell,ndim,dims);        
 
            delete[] dims_cell;
            delete[] ndim;
            delete[] dims[0];
            delete[] dims;

    } else if(type == 2 && par->sim_only_inc == 0){

        par->C_ast = new double*[par->T*par->Nd_max];
        par->M_ast = new double*[par->T*par->Nd_max];
        par->V_ast = new double*[par->T*par->Nd_max];        

        for(int t = 0; t < par->T;      t++){
        for(int d = 0; d < par->Nd_max; d++){          

            int i_cell = index::d2(t,d,par->Nd_max);
            par->C_ast[i_cell] = (double*) mxGetPr(mxGetCell(mxGetField(prhs[1],0,"C_ast"),i_cell));
            par->M_ast[i_cell] = (double*) mxGetPr(mxGetCell(mxGetField(prhs[1],0,"M_ast"),i_cell));
            par->V_ast[i_cell] = (double*) mxGetPr(mxGetCell(mxGetField(prhs[1],0,"V_ast"),i_cell));                         

        } }

    }

    ///////////////////////////
    // 4. outputs - simulate //
    ///////////////////////////

    if(type == 2){

        // a. struct
        const char *field_names[] = {"d", "M", "P", "Y", "N", "C", "A", "V", "MPC"};
        
        int num_fields = sizeof(field_names)/sizeof(*field_names);
        plhs[0] = mxCreateStructMatrix(1, 1, num_fields, field_names);
        auto sim_struct = plhs[0]; 

        // b. dimensions
        int ndim = 2;
        auto dims = new size_t[2];

            dims[0] = par->simN;
            dims[1] = par->simT;

        // c. elements
        par->d = mymex::set_field_int(sim_struct,"d",ndim,dims);
        par->M = mymex::set_field_double(sim_struct,"M",ndim,dims);
        par->N = mymex::set_field_double(sim_struct,"N",ndim,dims);
        par->P = mymex::set_field_double(sim_struct,"P",ndim,dims);
        par->Y = mymex::set_field_double(sim_struct,"Y",ndim,dims);        
        par->C = mymex::set_field_double(sim_struct,"C",ndim,dims);
        par->A = mymex::set_field_double(sim_struct,"A",ndim,dims);
        par->V = mymex::set_field_double(sim_struct,"V",ndim,dims);
        par->MPC = mymex::set_field_double(sim_struct,"MPC",ndim,dims);

            delete[] dims;

    }

        if(type == 1){
            logs::solve(0," output allocated\n");
        }
}


////////////////
// 3. destroy //
////////////////

void destroy(par_struct *par, int type){

    if(par->LRT == 1){
        delete[] par->i_d_plus;
        delete[] par->i_d_plus_p;
        delete[] par->i_d_plus_p_cumsum;
    }

    if(type == 2 && par->sim_only_inc == 0){
        delete[] par->C_ast;
        delete[] par->M_ast;
        delete[] par->V_ast;                
    }
}

} // namespace