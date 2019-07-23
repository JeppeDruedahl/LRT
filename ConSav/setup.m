classdef setup
methods(Static)

function [par] = LRT()
    
    par = struct();
    par.prefix = 'LRT';

    % a. grids
    par.NA = 150;
    par.NM = par.NA+200;    
    par.NN = 100;
    par.NP = 1;

        % grid specifications
        par.Amin = 1e-6;
        par.Amax = 2;
        par.phi_A = 1.2;
        
        par.Nmin = 0;
        par.Nmax = 2;
        par.phi_N = 1.1;
        
        par.Pmin = 0.1;
        par.Pmax = 3;
        par.phi_P = 1.1;
        
        par.Mmin = 1e-6;
        
    % b. demographics
    par.T = 31;
    par.age_min = 30;
    
    % c. preferences
    par.epstein_zin = 0;
    par.beta = 0.96;
    par.rho = 2.0;
    par.sigma = 2.0;
    par.zeta = 5.0;

    % d. old income process
    par.sigma_psi = 0.0;
    par.sigma_xi = 0.0;

        par.Npsi = 1;
        par.Nxi = 1;
    
    % f. assets
    par.R = 1.03;
    par.Rb = 1.04;
    par.kappa = 0.05;

    % g. technical
    par.LRT = 1; 
    par.tmin = 0;
    par.inceq = 0.0;
    par.sim_only_inc = 0;
        
    % h. simulate
    par.simT = par.T-1;
    par.simN = 500000;
    par.A_lag_ini = 0.0;
    par.P_ini_std = 0.0;

    % load data
    par = model.load_data(par);
    
end
function par = simple_LRT()

    par = setup.LRT();
    par.prefix = 'simple_LRT';
    par.A_lag_ini = 3.0;
    
        par.epstein_zin = 0;
        par.beta = 0.96;
        par.R = par.beta^-1;
        par.rho = 2.0;  
        par.sigma = 2.0;
        par.zeta = 1.0;
        par.kappa = 0.0;

    % a. number of elements 
    par.Nd = int32(10*ones(par.T,1));
    par.Nd_max = max(par.Nd);

    % b. income matrix
    par.grid_Y = nan(par.Nd_max,par.T);
    for t = 1:par.T
    for i = 1:par.Nd(t)
        if t == par.T
            par.grid_Y(i,t) = 0.0;
        else
            par.grid_Y(i,t) = 1.0;
        end
    end
    end
       
    % c. transition matrices
    par.Nd_plus = int32(zeros(par.Nd_max,par.T));
    par.i_d_plus = cell(par.Nd_max,par.T);
    par.i_d_plus_p = cell(par.Nd_max,par.T);
    par.i_d_plus_p_cumsum = cell(par.Nd_max,par.T);

    for t = 1:par.T-1
    for i = 1:par.Nd_max
    
        par.Nd_plus(i,t) = par.Nd(t+1);   
        par.i_d_plus{i,t} = int32(0:9); % remember c-indexing
        par.i_d_plus_p{i,t} = 0.1*ones(par.Nd(t+1),1);   
        par.i_d_plus_p_cumsum{i,t} = cumsum(par.i_d_plus_p{i,t});

    end    
    end

    % d. growth income growth rates
    par.Gamma = ones(par.T,1);

    % e. initial groups
    par.prob_G_ini = 0.1*ones(par.Nd(1),1);   

end
function par = PIH()

    par = setup.LRT();
    par.prefix = 'PIH';
    par.LRT = 0;
    par.A_lag_ini = 3.0;
    
        par.Nd_max = 1;
        par.Nd = int32(ones(par.T,1));
        par.NP = 100;

        par.epstein_zin = 0;
        par.beta = 0.96;
        par.R = par.beta^-1;
        par.rho = 2.0;   
        par.sigma = 2.0;
        par.zeta = 1.0;
        par.kappa = 0.0;

end
function par = noLRT()

    par = setup.LRT();
    par.LRT = 0;
    par.prefix = 'noLRT';

        par.Nd_max = 1;
        par.Nd = int32(ones(par.T,1));
        par.NP = 60;
        
    % PT income process
    tbl = readtable('../output/PT_estimates.txt','ReadVariableNames',false);
    par.sigma_psi = tbl.Var1(1);
    par.sigma_xi = tbl.Var1(2);
    par.P_ini_std  = tbl.Var1(3);

        par.Npsi = 8;
        par.Nxi = 8;

end
function par = no_risk()

    par = setup.LRT();
    par.LRT = 0;
    par.prefix = 'no_risk';

        par.Nd_max = 1;
        par.Nd = int32(ones(par.T,1));
        par.NP = 60;
        
    % PT income process
    par.P_ini_std = 0.0;   

    par.sigma_xi = 0.0;
    par.sigma_psi = 0.0;

        par.Npsi = 6;
        par.Nxi = 6;

end

end
end