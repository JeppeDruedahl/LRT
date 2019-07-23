models = {'LRT','LRT_ez_rho2','LRT_ez_rho3','LRT_ez_rho4'};

table = cell(1,1);
for i = 1:numel(models)
   
    % a. load
    model = models{i};
    load(sprintf('../output/ConSavEstimates/%s.mat',model))
    
    % b. preferences
    j = 1;
    table{j,i} = sprintf('%5.3f',par.sigma); j=j+1;
    table{j,i} = sprintf('%5.3f',par.rho); j=j+1;
    table{j,i} = ''; j=j+1;
    table{j,i} = sprintf('%5.3f',par.beta); j=j+1;
    table{j,i} = sprintf('%5.3f',par.zeta); j=j+1;
    table{j,i} = ''; j=j+1;

    % c. welfare cost
    table{j,i} = sprintf('%5.3f',par.inceq_est); j=j+1;  
    table{j,i} = ''; j=j+1;

    % d. fit
    if isfield(par,'preferences_fval')
        table{j,i} = sprintf('%5.3f',par.preferences_fval/1000); j=j+1;   
    end    

    
end
rownames = {'Intertemporal subsitution, $\sigma$','Risk aversion, $\rho$','',...
            'Discount factor, $\beta$','Post-retirement saving motive, $\zeta$','',...
            'Welfare cost of income risk, $\delta$','',...
            'Fit of $A_t$ (MSE)'};
funs.printtable(rownames,table,'figs_tabs/main_table.tex');