clc;
clear;
close all;
funs.layout();

%% 1. MEX

threads = 72;
model.mex_solve(threads);
model.mex_simulate(threads);
   
%% 2. settings

% a. which
do_simple = 0;
do_noLRT = 0; % 1
do_LRT = 1;
do_CRRA = 0; % 0
do_ez = 1; 
do_extra_noLRT = 0; % 1
rhos_noLRT = [];
rhos_LRT = [4]; %[2,3,4]; TEMP TEMP TEMP

% b. what
do_estimate = 1;
do_pol = 0;
do_lcp = 1;
find_inc_eq = 1;

%% 3. test I: PIH (should have flat consumption profile)

if do_simple 
    
    [par_PIH,sol_PIH,sim_PIH] = model.all('PIH',{},{},0,do_pol,do_lcp,find_inc_eq);

        % epstein-zin: should give the same because there is no uncertainty
        if do_lcp
            parnames = {'prefix','epstein_zin','rho'};
            parvals = {'PIH_ez',1,4};
            model.all('PIH',parnames,parvals,0,do_pol,do_lcp,find_inc_eq);    
        end
        
end

%% 4. test II: simple_LRT (should have flat consumption profile)

if do_simple
    
    [par_sLRT,sol_sLRT,sim_sLRT] = model.all('simple_LRT',{},{},0,do_pol,do_lcp,find_inc_eq);

        % epstein-zin: should give the same because there is no uncertainty
        if do_lcp
            parnames = {'prefix','epstein_zin','rho'};
            parvals = {'simple_LRT_ez',1,4};
            model.all('simple_LRT',parnames,parvals,0,do_pol,do_lcp,find_inc_eq);
        end
        
end
        
%% 5. noLRT model

if do_noLRT
    
    if do_CRRA
        [par_noLRT,sol_noLRT,sim_noLRT] = model.all('noLRT',{},{},do_estimate,do_pol,do_lcp,find_inc_eq);
    end
    
    % epstein-zin: should give the same because rho = sigma
    if do_ez 
        
        if do_CRRA
            parnames = {'prefix','epstein_zin','beta','kappa','zeta'};
            parvals = {'noLRT_ez',1,par_noLRT.beta,par_noLRT.kappa,par_noLRT.zeta};
            if do_estimate == 2
                model.all('noLRT',parnames,parvals,2,do_pol,do_lcp,find_inc_eq);
            else
                model.all('noLRT',parnames,parvals,0,do_pol,do_lcp,find_inc_eq);
            end
        end
        
        for rho = rhos_noLRT
            parnames = {'prefix','epstein_zin','sigma','rho'};
            parvals = {sprintf('noLRT_ez_rho%d',rho),1,2/3,rho};
            model.all('noLRT',parnames,parvals,do_estimate,do_pol,do_lcp,find_inc_eq);
        end
        
    end
    
    % extra
    if do_extra_noLRT == 1
        parnames = {'prefix','sigma_psi','sigma_xi','P_ini_std'};
        parvals = {'noLRT_low_vars',0.08, 0.08, 0.20};        
        model.all('noLRT',parnames,parvals,do_estimate,do_pol,do_lcp,find_inc_eq);
    end
    
end

%% 6. LRT model

if do_LRT
    
    if do_CRRA
        [par_LRT,sol_LRT,sim_LRT] = model.all('LRT',{},{},do_estimate,do_pol,do_lcp,find_inc_eq);
    end
    
        % epstein-zin: should give the same because rho = sigma
        if do_ez
            
            if do_CRRA
                parnames = {'prefix','epstein_zin','beta','kappa','zeta'};
                parvals = {'LRT_ez',1,par_LRT.beta,par_LRT.kappa,par_LRT.zeta};
                if do_estimate == 2
                    model.all('LRT',parnames,parvals,2,do_pol,do_lcp,find_inc_eq);
                else
                    model.all('LRT',parnames,parvals,0,do_pol,do_lcp,find_inc_eq);
                end                
            end
            
            for rho = rhos_LRT
                parnames = {'prefix','epstein_zin','sigma','rho'};
                parvals = {sprintf('LRT_ez_rho%d',rho),1,2/3,rho};
                model.all('LRT',parnames,parvals,do_estimate,do_pol,do_lcp,find_inc_eq);
           end
        end
        
end

%% 7. move figures and tables

write_table
move_figs_tabs

%% 8. clean up

rmdir('figs_tabs')
delete('log*.txt')
delete('*.mexw64')