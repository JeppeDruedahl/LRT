namespace utilities {

////////////////
// 1. utility //
////////////////

double value_of_choice(double C, double W, par_struct *par)
{

    if(par->epstein_zin == 0){
        double u = pow(C,1.0-par->rho)/(1.0-par->rho);
        return (1.0-par->beta)*u +par->beta*W;
    } else {
        double u1 = (1.0-par->beta)*pow(C,1.0-par->sigma);
        double u2 = par->beta*pow(W,1.0-par->sigma);
        return pow(u1+u2,1.0/(1.0-par->sigma));
    }

}

double marg_u(double C, par_struct *par)
{
    if(par->epstein_zin == 0){    
        return pow(C,-par->rho);
    } else {
        return pow(C,-par->sigma);        
    }
}

double inv_marg_u(double u, par_struct *par)
{
    if(par->epstein_zin == 0){      
        return pow(u,-1.0/par->rho);
    } else {
        return pow(u,-1.0/par->sigma);
    }
}

}