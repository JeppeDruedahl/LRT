classdef model
methods(Static)

function mex_solve(threads)
       
    str = sprintf('mex cfuncs/solve_model.cpp -DMAXTHREADS=%d',threads);
    str = sprintf('%s %s',str,' COMPFLAGS="$COMPFLAGS /Ox /Wall /openmp"');

    eval(str);
    
end
function mex_simulate(threads)
    
    str = sprintf('mex cfuncs/simulate_model.cpp -DMAXTHREADS=%d',threads);
    str = sprintf('%s %s',str,' COMPFLAGS="$COMPFLAGS /Ox /Wall /openmp"');

    eval(str);
    
end
function par = load_data(par)
    
    % a. number of elements 
    tbl = readtable('../data/num_leafs.txt','ReadVariableNames',false);
    par.Nd = int32([tbl.Var1;1]);
    par.Nd_max = max(par.Nd);

    % b. income matrix
    par.grid_Y = nan(par.Nd_max,par.T);
    for t = 1:par.T        
        if t == par.T
            par.grid_Y(1,t) = 0.0;
        else
            tbl = readtable(sprintf('../data/ypred_G_t%d.txt',t-1),'ReadVariableNames',false);        
            par.grid_Y(1:par.Nd(t),t) = exp(tbl.Var1);      
        end
    end
       
    % c. transition matrices
    par.Nd_plus = int32(zeros(par.Nd_max,par.T));
    par.i_d_plus = cell(par.Nd_max,par.T);
    par.i_d_plus_p = cell(par.Nd_max,par.T);
    par.i_d_plus_p_cumsum = cell(par.Nd_max,par.T);

    t = par.T-1;
    for i = 1:par.Nd(t)
        par.Nd_plus(i,t) = 1;             
        par.i_d_plus{i,t} = int32(0);
        par.i_d_plus_p{i,t} = 1.0;   
        par.i_d_plus_p_cumsum{i,t} = 1.0;
    end

    for t = 1:par.T-2

        tbl = readtable(sprintf('../data/trans_t%d.txt',t));
        trans = table2array(tbl);
        next = (1:par.Nd(t+1))-1; % remember c indexing
        
            % check dimensions
            assert(size(trans,1)==par.Nd(t));
            assert(size(trans,2)==par.Nd(t+1));
        
        for i = 1:par.Nd(t)        
            
            % i. transition vector
            I = isnan(trans(i,:)) == 0;
            trans_vec = trans(i,I);
            
            % ii. number of elements
            par.Nd_plus(i,t) = numel(trans_vec); 
            
            % iii. groups
            par.i_d_plus{i,t} = int32(next(I));
            
            % iv. probabilities
            par.i_d_plus_p{i,t} = trans_vec/sum(trans_vec);   
            par.i_d_plus_p_cumsum{i,t} = cumsum(par.i_d_plus_p{i,t});
            
        end
        
    end
        
    % e. initial groups
    tbl = readtable('../data/prob_G_ini.txt','ReadVariableNames',false);
    par.prob_G_ini = tbl.Var1;
    
    % f. mean and var of Y and logY
    tbl = readtable('../data/mean_Y.txt','ReadVariableNames',false);
    par.mean_Y = tbl.Var1;  

    tbl = readtable('../data/mean_Y_lev.txt','ReadVariableNames',false);
    par.mean_Y_lev = tbl.Var1; 
    
    tbl = readtable('../data/mean_logY.txt','ReadVariableNames',false);
    par.mean_logY = tbl.Var1;  
  
    tbl = readtable('../data/var_logY.txt','ReadVariableNames',false);
    par.var_logY = tbl.Var1;    
    
    % g. income growth rates
    par.Gamma = ones(par.T,1);
    for t = 1:(par.T-2)
        par.Gamma(t) = par.mean_Y(t+1)/par.mean_Y(t);
    end
    
    % h. retirement asssets
    par.p50_N_retire = 1212526*(1-0.4)/(1000*par.mean_Y_lev(1)); % KORA (2016, tabel 2.2)
    par.mean_N_retire = 1803548*(1-0.4)/(1000*par.mean_Y_lev(1)); % KORA (2016, tabel 2.2)

    % i. net worth
    tbl = readtable('../data/a_p50_p100.csv','ReadVariableNames',false);
    par.p50_A = tbl.Var1;
    par.p50_A = par.p50_A /par.mean_Y_lev(1); 

    tbl = readtable('../data/a_mean_p100.csv','ReadVariableNames',false);
    par.mean_A = tbl.Var1;
    par.mean_A = par.mean_A /par.mean_Y_lev(1); 
    
end
function par = create_grids(par)

    assert( (par.LRT == 1 && par.NP == 1) || par.LRT == 0 );
        
    % 1. grids
    par.grid_A = nan(par.NA,par.T);
    par.grid_P = nan(par.NP,par.T);
    par.grid_N = nan(par.NN,par.T);    
    for t = 1:par.T

        scale_P = prod(par.Gamma(1:t))*1.01^30;
        par.grid_P(:,t) = funs.nonlinspace(par.Pmin,par.Pmax,par.NP,par.phi_P)*scale_P;        
        
        scale_A = scale_P*(par.R+0.01)^30; 
        par.grid_A(:,t) = funs.nonlinspace(par.Amin,par.Amax,par.NA,par.phi_A)*scale_A;

        scale_N = scale_P*(par.Rb+0.01)^30;
        par.grid_N(:,t) = funs.nonlinspace(par.Nmin,par.Nmax,par.NN,par.phi_N)*scale_N;
    
    end
        
    % 2. shocks
    [psi, psi_w] = funs.GaussHermite_lognorm(par.sigma_psi,par.Npsi,1);
    [xi, xi_w] = funs.GaussHermite_lognorm(par.sigma_xi,par.Nxi,1);
        
    [psi, xi] = ndgrid(psi,xi);
    par.psi = psi(:);
    par.xi = xi(:);
    
    [psi_w, xi_w] = ndgrid(psi_w,xi_w);
    weights = psi_w.*xi_w;
    par.weights = weights(:);
            
    par.Nshocks_max = numel(par.weights);
            
        % convert into matrices
        par.psi = [repmat(par.psi,[1 par.T-1]) zeros(par.Nshocks_max,1)];
        par.xi = [repmat(par.xi,[1 par.T-1]) zeros(par.Nshocks_max,1)];
        par.weights = [repmat(par.weights,[1 par.T-1]) ones(par.Nshocks_max,1)];

        par.Nshocks =  int32([repmat(par.Nshocks_max,[1 par.T-1]) 1]);   

end
function sol = solve(par)

    sol = solve_model(par);

end
function par = setup_simulate(par)

    if par.LRT == 0
        par.d_ini = int32(zeros(par.simN,1));
    else
        d_ini = randsample(par.Nd(1),par.simN,'true',par.prob_G_ini);
        par.d_ini = int32(d_ini-1);
    end
    par.P_ini = exp(par.P_ini_std*randn(par.simN,1) - 0.5*par.P_ini_std^2);
    par.Y_ini = nan*ones(par.simN,1);    
    par.M_ini = par.R*par.A_lag_ini*ones(par.simN,1);
    par.N_ini = nan*par.Y_ini;

    par.p_sim = rand(par.simN,par.simT);

    par.psi_sim = exp(par.sigma_psi*randn(par.simN,par.simT) - 0.5*par.sigma_psi^2);
    par.xi_sim = exp(par.sigma_xi*randn(par.simN,par.simT) - 0.5*par.sigma_xi^2);
    
end
function [par] = updatepar(par,parnames,parvals)

    for i = 1:numel(parnames)           
        parname = parnames{i};
        if iscell(parvals)
            parval  = parvals{i};            
        else
            parval  = parvals(i);          
        end
        par.(parname) = parval;            
    end      

end 
function sim = simulate(par,sol)

    sim = simulate_model(par,sol);
    sim.logY = log(sim.Y);
    
end
function f = obj_kappa(x,par)
    
    par.kappa = x;

    par.sim_only_inc = 1;
    sim = model.simulate(par,nan);
    par.sim_only_inc = 0;

    %f = median(sim.N(:,end)) - par.p50_N_retire;
    f = mean(sim.N(:,end)) - par.mean_N_retire;

end
function par = estimate_kappa(par)

    goal = @(x) model.obj_kappa(x,par);
    options = optimoptions('fsolve','Display','none');
    [x,fval,exitflag,output] = fsolve(goal,0.1,options);

    par.kappa = x;
    par.kappa_fval = fval;

end
function f = obj_preferences(x,par)
    
    fprintf(' beta = %10.8f, zeta = %10.8f',x);    
    par.beta = x(1);
    par.zeta = x(2);
    sol = model.solve(par);
    sim = model.simulate(par,sol);
    
    f = 0;
    n = 0;
    for t = 1:par.T-1
        if isnan(par.p50_A(t)) == 0
            %f = f + 1000*(median(sim.A(:,t)) - par.p50_A(t)).^2;
            f = f + 1000*(mean(sim.A(:,t)) - par.mean_A(t)).^2;
            n = n + 1;
        end
    end
    f = f/n;

    fprintf(', f = %14.8f\n',f);

end
function par = estimate_preferences(par)

    goal = @(x) model.obj_preferences(x,par);

    options = optimoptions(@fmincon,'Display','none',...
        'StepTolerance',1e-4,'OptimalityTolerance',1e-4);
    lb = [0.85 0.50];
    ub = [0.99 10.00];
    x0 = [par.beta,par.zeta];
    [x,fval,exitflag,output] = funs.fmincon(goal,x0,lb,ub,options);

    par.beta = x(1);
    par.zeta = x(2);    
    par.preferences_fval = fval;
    fprintf('\n')

end
function V = calc_V(par,sim)

    if par.epstein_zin == 0
        V = mean(sim.V(:,1));
    else
        V = mean(sim.V(:,1).^(1-par.rho)).^(1/(1-par.rho));
    end

end
function f = obj_inceq(x,par,V_base)
    
    par.inceq = x;
    sol = model.solve(par);
    sim = model.simulate(par,sol);
    V = model.calc_V(par,sim);
   
    f = V-V_base;
    
end
function [par] = income_equivalent(par,sol,sim,do_lcp)

    % 1. baseline V
    V_base = model.calc_V(par,sim);

    % 1. no risk
    par_norisk = par;
    if par.LRT == 0
        
        par_norisk.sigma_psi = 0.0;
        par_norisk.sigma_xi = 0.0;
        par_norisk.P_ini_std = 0.0;
        
            par_norisk.Npsi = 1;
            par_norisk.Nxi = 1;

    else
        
        par_norisk.LRT = 0;
        
        par_norisk.Nd_max = 1;
        par_norisk.Nd = int32(ones(par.T,1));
        par_norisk.NP = 100;
        
        par_norisk.sigma_psi = 0.0;
        par_norisk.sigma_xi = 0.0;        
        par_norisk.P_ini_std = 0.0;

            par_norisk.Npsi = 1;
            par_norisk.Nxi = 1;
            
        par_norisk.Gamma = ones(par.T,1);
        for t = 1:(par.T-2)
            par_norisk.Gamma(t) = mean(sim.Y(:,t+1))/mean(sim.Y(:,t));
        end

    end

    par_norisk = model.create_grids(par_norisk);    
    par_norisk = model.setup_simulate(par_norisk);
   
    % 2. comparision
    if do_lcp
        sol_norisk = model.solve(par_norisk);        
        sim_norisk = model.simulate(par_norisk,sol_norisk);
        labels = {'full risk', 'no risk'};
        figs.lcp_compare('C','mean','$C_t$',par,sim,sim_norisk,'norisk',labels);
        figs.lcp_compare('A','mean','$A_t$',par,sim,sim_norisk,'norisk',labels);
        figs.lcp_compare('Y','mean','$Y_t$',par,sim,sim_norisk,'norisk',labels);
        figs.lcp_compare('logY','mean','$\log(Y_t)$',par,sim,sim_norisk,'norisk',labels);
        figs.lcp_compare('logY','var','var($\log(Y_t)$)',par,sim,sim_norisk,'norisk',labels);
    end
    
    % 3. find income equivalent
    goal = @(x) model.obj_inceq(x,par_norisk,V_base);
    options = optimoptions(@fsolve,'Display','none');
    [x,fval,exitflag,output] = fsolve(goal,log(0.9),options);
    
    par.inceq_est = 1.0-exp(x);
    par.inceq_fval = fval;
        
end
function [par,sol,sim] = all(par_str,parnames,parvals,do_estimate,do_pol,do_lcp,find_inc_eq)
        
    rng(1986)
    
    % 1. setup
    par = setup.(par_str);
    par = model.updatepar(par,parnames,parvals);    
    par = model.create_grids(par);
    par = model.setup_simulate(par);

         % clean up
        folder = ['figs_tabs\' par.prefix];
        if exist(folder,'dir') > 0
            delete(sprintf('%s/*.png',folder));
            delete(sprintf('%s/*.pdf',folder));            
            delete(sprintf('%s/*.txt',folder));        
        else
            mkdir(folder);    
        end
        
        % print
        fprintf('\n\n');
        fprintf('**********************\n');
        fprintf('** %-16s **\n',par.prefix);
        fprintf('**********************\n\n');
         
    % 2. estimate
    if do_estimate == 1

        % a. kappa
        t0 = tic;
        par = model.estimate_kappa(par);
        estimate_time = toc(t0);
        fprintf(' estimate time: %10.2f secs\n',estimate_time)
        fprintf('kappa estimate: %10.4f [fval = %5.4f]\n\n',par.kappa,par.kappa_fval);            

        % b. beta
        t0 = tic;
        par = model.estimate_preferences(par);
        estimate_time = toc(t0);
        fprintf(' estimate time: %10.2f secs\n',estimate_time)
        fprintf(' beta estimate: %10.4f\n',par.beta); 
        fprintf(' zeta estimate: %10.4f\n',par.zeta);         
        fprintf('          fval: %10.4f\n\n',par.preferences_fval);
        
    elseif do_estimate == 2
        
        temp = par;
        load(sprintf('../output/ConSavEstimates/%s',par.prefix))
        
        temp.kappa = par.kappa;
        if isfield(par,'kappa_fval')
            temp.kappa_fval = par.kappa_fval;
        end
        
        temp.beta = par.beta;
        temp.zeta = par.zeta;
        if isfield(par,'preferences_fval')
            temp.preferences_fval = par.preferences_fval;
        end
        
        temp.inceq_est = par.inceq_est;
        temp.inceq_fval = par.inceq_fval;

        par = temp;
        
    end

    % 3. solve  
    t0 = tic;
    sol = model.solve(par);
    sol_time = toc(t0);
    fprintf('solve time:    %10.2f secs\n',sol_time)
    
        if do_pol
            figs.all_polfuncs(par,sol);
        end
        
    % 4. simulate
    t0 = tic;
    sim = model.simulate(par,sol);
    sim_time = toc(t0);
    fprintf('simulate time: %10.2f secs\n',sim_time)
       
        if do_lcp
            
            figs.lcp('C','mean','$C_t$',par,sim);         
            figs.lcp('A','mean','$A_t$',par,sim); 
            figs.lcp('N','mean','$N_t$',par,sim);        
            figs.lcp('MPC','mean','MPC',par,sim);    

            labels = {'model','data'};
            figs.lcp_compare('Y','mean','$Y_t$',par,sim,par.mean_Y,'data',labels);
            figs.lcp_compare('logY','mean','$log(Y_t)$',par,sim,par.mean_logY,'data',labels);
            figs.lcp_compare('logY','var','$var(log(Y_t))$',par,sim,par.var_logY,'data',labels);  
            figs.lcp_compare('A','mean','$A_t$',par,sim,par.mean_A,'data',labels);
            
            for perc = [5 95]
                figs.lcp('Y',perc,'',par,sim);
                figs.lcp('A',perc,'',par,sim);                
                figs.lcp('N',perc,'',par,sim);                                
            end
            
            if par.Nd_max > 1
                figs.d_hist([1 2 10 20 30],par,sim);
            end
            
        end
        
    % 5. income equivalent
    if find_inc_eq == 1 && do_estimate ~= 2
        t0 = tic;
        par = model.income_equivalent(par,sol,sim,do_lcp);
        inceq_time = toc(t0);
        fprintf('inceq time:    %10.2f secs\n\n',inceq_time);
        fprintf('income equivalent: %10.4f [fval = %6.4f]\n\n',par.inceq_est,par.inceq_fval);
    else
        fprintf('\n')
    end
    
    % 6. print some info
    vars = {'C','A'};
    for t = [20]
    for i = 1:numel(vars)
        fprintf('%3s in t = %2d:  %10.4f (mean) %10.4f (var)\n',...
            vars{i},t,mean(sim.(vars{i})(:,t)),var(sim.(vars{i})(:,t)));
    end
    end
    fprintf('\n')
    
    % 7. save
    par.MPC_mean = mean(funs.vec(sim.MPC));
    par.MPC_p50 = median(funs.vec(sim.MPC));    
    
    t0 = tic;
    par.grid_Y = [];
    par.grid_N = [];
    par.grid_A = [];
    par.grid_P = [];
    par.i_d_plus = [];
    par.i_d_plus_p = [];        
    par.i_d_plus_p_cumsum = [];
    par.d_ini = [];
    par.P_ini = [];
    par.Y_ini = [];       
    par.M_ini = [];
    par.N_ini = [];                
    par.p_sim = [];
    par.psi_sim = [];
    par.xi_sim = [];
    
    folder = '../output/ConSavEstimates' ;
    if exist(folder,'dir') == 0
        mkdir(folder);
    end

    save(sprintf('../output/ConSavEstimates/%s.mat',par.prefix),'par','-v7.3');
    save_time = toc(t0);
    fprintf(' save time:    %10.2f secs\n\n',save_time);

end

end
end