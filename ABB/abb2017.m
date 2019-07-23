% Class function for Arellano, Blundell, and Bonhomme (2017; ECTA)
%
% Code works with two basic structs: data and par. 
%   data: contains the raw income data as well as convenient
%   transformations needed in the various estimation functions. 
%   par: model parameters and technical flags. 
% 
% Parameter estimates are saved in three structs: par.eta_0, par.eta_t, and
% par.eps. These are objects of the type "abb_par" (see abb_par.m). 
%
%
% Parameters: there are three equations corresponding to eta_0, eta_t and
% eps. Each have an associated quantile regression with an x and a y
% variable, and parameters for each segment of [0;1] (indexed by jtau)  
% - par.eta_0:  x = MatAGE1_tot (initial age), y = Matdraw(:,1) 
% - par.eta_t:  x = Matdraw_lag,               y = Matdraw_t 
% - par.eps:    x = MatAGE_tot,                y = Ytot_t-Matdraw_tot 
%

classdef abb2017
    
    methods (Static)
                
        function [data, par] = initialize(Y, AGE, T, K1, K2, K3, K4, KPolyAgeResid, Vectau, var_prop, draws, maxiter)
            % initialize: form the data and par structs based on raw income
            % data and a range of settings. 
            %
            % INPUTS: 
            %   Y: N*T_ matrix of 
            %   T: how many years to use for each obs. (will be randomly
            %      drawn). T <= T_ must hold 
            %   Vectau: Vector over [0;1] where the spline grid points will
            %      be 
            %   draws: (integer) 
            %   maxiter: (integer) 
            
            [N,T_] = size(Y);
            assert( all(size(AGE) == size(Y)), 'The matrices Y and AGE must have same dimensions'); 
            
            % subset
            if T<T_
                fprintf('randomly selecting a subset of ages for each individual\n');
                Y_full = Y;
                AGE_full = AGE;
                Y = nan(N,T);
                AGE = nan(N,T);
                ages = randi(T_-T, N, 1); % if ages start too late, they end up out of bounds
                for i=1:N
                    i_row =  ages(i):ages(i)+T-1;
                    Y(i,:) = Y_full(i,i_row) ;
                    AGE(i,:) = AGE_full(i,i_row);
                end
            elseif T == T_ 
                fprintf('Keeping the full dataset\n');
            else 
                error('Asking for T=%d years of income data but only %d available in input data.', T, T_)
            end
            
            assert(all(size(Y) == [N,T]), 'Unexpected dimension of Y: subsampling above must have failed');
            
            meanAGE=mean(AGE(:));
            stdAGE=std(AGE(:));
            
            % --- residuzlize Y against age ---
            % Regression of earnings on a fourth order Hermite polynomial in age
            XX = abb2017.hermite_basis_univariate((AGE(:)-meanAGE)/stdAGE, KPolyAgeResid-1);
            coeff=pinv(XX)*Y(:); % equivalent of (XX'*XX)\XX'*Y(:)
            Ypred=XX*coeff;
            Y = Y - reshape(Ypred,N,T);
            
            % show results
            Vectage=(25:1:60)';
            XX1 = abb2017.hermite_basis_univariate((Vectage-meanAGE)/stdAGE, KPolyAgeResid-1);
            
            % uncomment to show a plot of residualized income over age s
            %plot(Vectage,XX1*coeff); xlabel('Age (residualized)'); ylabel('Predicted income');
            
            meanY=mean(Y(:));
            stdY=std(Y(:));
            
            % reshape 
            Y_t   = reshape(Y(:,2:T),    N*(T-1), 1);
            AGE_t = reshape(AGE(:, 2:T), N*(T-1), 1);
            Ylag  = reshape(Y(:,1:T-1),  N*(T-1), 1);
            
            % regressors for qreg of Y_t on Y_t-1 and age (polynomials in each)
            % dimension = (N*(T-1), (K1+1)*(K2+1));
            MatYlag = abb2017.hermite_basis_bivariate((Ylag-meanY)/stdY, (AGE_t-meanAGE)/stdAGE, K1, K2);
            
            Mdraws = 1; % from the original code
            Ytot_t=repmat(Y(:,1), Mdraws);
            for j=2:T
                Ytot_t=[Ytot_t;repmat(Y(:,j), Mdraws)];
            end
            
            % when Mdraws == 1, this is == AGE_t
            MatAGE_t_tot=repmat(AGE_t(1:N,:), Mdraws, 1);
            for j=2:T-1
                MatAGE_t_tot=[MatAGE_t_tot;repmat(AGE_t((j-1)*N+1:j*N,:), Mdraws)];
            end
            
            MatAGE1  = abb2017.hermite_basis_univariate((AGE(:,1)-meanAGE)/stdAGE, K3); % dim = N by K3+1
            MatAGE_t = abb2017.hermite_basis_univariate((AGE(:)-meanAGE)/stdAGE, K4);   % dim = N*T by K4+1
            
            MatAGE1_tot=repmat(MatAGE1, Mdraws);
            MatAGE_tot=repmat(MatAGE_t(1:N,:), Mdraws, 1);
            for t=2:T
                rows_t = (t-1)*N+1:t*N;
                MatAGE_tot=[MatAGE_tot; ...
                    repmat(MatAGE_t(rows_t,:), Mdraws, 1)];
            end
            
            data.Y = Y;
            data.AGE = AGE;
            data.MatAGE_t = MatAGE_t;
            data.MatYlag = MatYlag;
            data.Y_t = Y_t;
            % and the Mdraws repeats
            data.MatAGE_t_tot = data.MatAGE_t; % assuming that Mdraws == 1!!
            data.Ytot_t = data.Y;
            data.MatAGE1 = MatAGE1;
            data.MatAGE1_tot = MatAGE1; % duplicate
            
            par.Vectau = Vectau;
            par.meanAGE = meanAGE;
            par.stdAGE = stdAGE;
            par.meanY = meanY;
            par.stdY = stdY;
            par.K1 = K1;
            par.K2 = K2;
            par.K3 = K3;
            par.K4 = K4;
            par.KPolyAgeResid = KPolyAgeResid; 
            par.maxiter = maxiter;
            par.coeff_Y_on_age_standardized = coeff;
            
            assert(~isempty(K3), 'Cannot find global K3')
            
            % parameters
            par.N = N;
            par.T = T;
            
            par = abb2017.update_par_mcmc_from_globals(par);
            data = abb2017.precompute(data,par);
            par = abb2017.precompute_par(par);
            par.var_prop = var_prop;
            par.draws = draws;
            
            % initialize chain of estimates with empty structs
            par.estimates_chain.eta_0(draws,1) = abb_par();
            par.estimates_chain.eta_t(draws,1) = abb_par();
            par.estimates_chain.eps(draws,1) = abb_par();
            
            par.USEFMIN = false; % not implemented
            %if par.USEFMIN
            %    par.initial = struct();
            %    par.initial.eta_0 = abb_par(Resqinit_e0, b1init_e0, bLinit_e0);
            %    par.initial.eta_t = abb_par(Resqinit, b1init, bLinit);
            %    par.initial.eps   = abb_par(Resqinit_eps, b1init_eps, bLinit_eps);
            %end
            
        end
        
        function par = estimate(data, par, USEWAITBAR)
            % Estimate the model with MCMC
            % 
            % OUTPUT: 
            %   Populates par.estimates_chain, which contains arrays of
            %   each of the three structs for eta_t, eta_0 and eps (which
            %   are of the class "abb_par". 
            %
            
            maxiter = par.maxiter; 
            
            init = randn(par.N,par.T);
            Obj_chain = [abb2017.density(init, data, par), zeros(par.N,par.draws-1)];
            
            % initialize
            myT = nan(maxiter, 1);
            Nu_chain = repmat({init}, par.draws,1);
            acceptrate = zeros(par.T,par.draws);
            mat_b=zeros(maxiter,6);
            mat_lik=zeros(maxiter,1);
            
            
            if USEWAITBAR
                h = waitbar(0, 'First iteration');
            end
            
            for iter=1:maxiter
                tic
                
                % 0. initial
                if iter>1
                    % initial values: from last iteration
                    Obj_chain(:,1) = Obj_chain(:,end);
                    Nu_chain{1} = Nu_chain{par.draws};
                end
                
                % 1. E-step
                for j = 2:par.draws
                    [Nu_chain{j}, Obj_chain(:,j), acceptrate(:,j)] = abb2017.update_chain(data, par, Nu_chain{j-1}, Obj_chain(:,j-1));
                end
                
                % 2. M-step
                [par, Matdraw] = abb2017.m_step(Nu_chain{par.draws}, data, par);
                
                % 3. update estimates chain
                par.estimates_chain = abb2017.update_estimates_chain(par, iter);
                
                % 4. print iteration info
                mat_lik(iter) = mean(log(abb2017.density(Matdraw, data,par)));
                abb2017.print_avg_acceptrate(acceptrate);
                fprintf('Full likelihood (iter %d): %8.4f.\n', iter, mat_lik(iter));
                % abb2017.quick_persistence_computation_avg(Nu_chain{end}, data, par);
                % abb2017.print_parameters(par);
                
                % 5. time to completion info
                myT(iter) = toc;
                avg_sec_per_iter = mean(myT(1:iter));
                now = datetime('now');
                E_finish = now + seconds(avg_sec_per_iter*(maxiter-iter));
                E_finish.Format = 'dd-MMM HH:mm';
                msg = sprintf('%d/%d: %5.1fs (avg. %5.1fs/it). ETA: %s.\n', iter, maxiter, myT(iter), avg_sec_per_iter, E_finish);
                if USEWAITBAR
                    waitbar(iter/maxiter, h, msg);
                else
                    fprintf(msg);
                end
                
            end
            
            if USEWAITBAR
                close(h);
            end
            
            
        end
        
        function [par] = estimate_initial_parameters(data, par)
            % Estimates initial parameters
            % RETURNS: populates the fields: 
            %   par.eta_t 
            %   par.eta_0 
            %   par.eps 
            
            fprintf('Estimating initial values.\n');
            
            Ntau = numel(par.Vectau); 
            
            % 1. eta_t equation. QREG of Y_t+N(0,1) on Y_t-1
            Resqinit=zeros((par.K1+1)*(par.K2+1),Ntau);
            fprintf('eta_t... ');
            for jtau=1:Ntau
                tau=par.Vectau(jtau);
                beta1=rq(data.MatYlag,data.Y_t+randn(par.N*(par.T-1),1),tau);
                Resqinit(1:(par.K1+1)*(par.K2+1),jtau)=beta1;
            end
            
            % 2. eta_0: initial condition given initial age
            
            % 2.a check that there is variation in initial age. Otherwise,
            % we cannot estimate a regression for eta_0 on a polynomial in
            % age. 
            assert(numel(unique(data.AGE(:,1)))>1, 'Without variation in initial age, this will not work');
            
            % 2.b estimate 
            Resqinit_e0=zeros(par.K3+1,Ntau);
            fprintf('eta_0... ');
            for jtau=1:Ntau
                tau=par.Vectau(jtau);
                beta1=rq(data.MatAGE1,data.Y(:,1)+randn(par.N,1),tau);
                Resqinit_e0(1:par.K3+1,jtau)=beta1;
            end
            
            % 3. eps: density of the error term
            Resqinit_eps=zeros(par.K4+1,Ntau);
            fprintf('eps... ');
            for jtau=1:Ntau
                tau=par.Vectau(jtau);
                beta1=rq(data.MatAGE_t,data.Y(:)+randn(par.N*par.T,1),tau);
                Resqinit_eps(1:par.K4+1,jtau)=beta1;
            end            
            
            fprintf('done!\n');
            
            % 4. initial Laplace parameter values
            b1      = 10;
            bL      = 10;
            b1_e0   = 10;
            bL_e0   = 10;
            b1_eps  = 10;
            bL_eps  = 10;

            % 5. populate par struct 
            par.eta_t.Param = Resqinit;
            par.eta_0.Param = Resqinit_e0; 
            par.eps.Param   = Resqinit_eps; 
            par.eta_t.b1    = b1;
            par.eta_t.bL    = bL;
            par.eta_0.b1    = b1_e0;
            par.eta_0.bL    = bL_e0;
            par.eps.b1      = b1_eps;
            par.eps.bL      = bL_eps;
            
        end
        
        function [data,par] = repackage_globals()
            % legacy function to convert from globals to data-par format
            
            global Y N T Ntau Vectau Resqinit_eps Resqinit Resqinit_e0 b1_eps bL_eps b1 bL b1_e0 bL_e0 K1 K2 K3 K4 AGE MatAGE_t MatAGE1  meanAGE stdAGE meanY stdY
            
            data = struct();
            
            data.Y = Y;
            data.AGE = AGE;
            data.MatAGE_t = MatAGE_t;
            data.MatAGE_t_tot = data.MatAGE_t; % assuming that Mdraws == 0!!
            data.Ytot_t = data.Y;
            data.MatAGE1_tot = MatAGE1;
            
            par.Vectau  = Vectau;
            par.meanAGE = meanAGE;
            par.stdAGE  = stdAGE;
            par.meanY   = meanY;
            par.stdY    = stdY;
            
            par.K1 = K1;
            par.K2 = K2;
            par.K3 = K3;
            par.K4 = K4;
            
            assert(~isempty(K3), 'Cannot find global K3')
            
            % parameters
            par.N = N;
            par.T = T;
            
        end
        
        function Ysim = simulate_data(Nsim, ages_sim, data, par)
            % simulate_data: simulates a balanced panel 
            % 
            % INPUTS: 
            %   Nsim: integer, # of individuals to simulate
            %   ages_sim: par.T-vector of corresponding ages 
            %
            % OUTPUT: 
            %   Ysim: Nsim*par.T matrix of simulated income histories. 
            %            
            
            assert( isscalar(Nsim), 'Nsim must be a scalar'); 
            assert( (size(ages_sim,1)==1) & (size(ages_sim,2)>1), 'ages_sim must be a col vector.'); 
            assert( numel(ages_sim) == par.T , sprintf('Must simulate a full T-panel (%d ages).', par.T));
            
            % prepare
            data2 = data;
            par2  = par;
            par2.N = Nsim;
            data2.AGE = repmat(ages_sim, par2.N, 1);
            par2.T = size(data2.AGE,2);
            data2 = abb2017.precompute(data2, par2);
            par2 = abb2017.precompute_par(par2);
            
            % simulate data
            rng('default');
            pmin_vec = .3*ones(par2.T,1)*min(data.Y(:));
            pmax_vec = 3*ones(par2.T,1)*max(data.Y(:));
            data2.V_draw_eta = rand([par2.N, par2.T]);
            data2.V_draw_eps = rand([par2.N, par2.T]);
            [Ysim2, eta2, eps2] = abb2017.simulate(data2, par2, pmin_vec, pmax_vec);
            
            % add back age effects (so that means will fit)
            XX = abb2017.hermite_basis_univariate( (data2.AGE(:)-par.meanAGE)/par.stdAGE , par.KPolyAgeResid-1);
            Ysim = Ysim2(:) + XX*par.coeff_Y_on_age_standardized;
            eta3 = eta2(:) + XX*par.coeff_Y_on_age_standardized;
            Ysim = reshape(Ysim, size(Ysim2));
            eta3 = reshape(eta3, size(eta2));
        end
        
        function simulate_and_compare_AR1_regression(data, par)
            % simulate dataset with same dimensions
            pmin_vec = [nan; 3*min(data.Y(:,2:end))']; % why min*3?
            pmax_vec = [nan; 3*max(data.Y(:,2:end))'];
            rng('default');
            data.V_draw_eta = unifrnd(0,1,par.N,par.T);
            data.V_draw_eps = unifrnd(0,1,par.N,par.T);
            [Ysim, eta, eps] = abb2017.simulate(data,par,pmin_vec,pmax_vec);
            
            fprintf('--- AR(1) regressions --- \n');
            fprintf('Below are shown coefficients from a simple AR(1) regression \nfirst for the simulated and then for the real data.\n'); 
            
            % ar(1) regs
            fprintf('Simulated data\n');
            Yt = Ysim(:,2:end);
            Yt_1 = Ysim(:,1:end-1);
            regress(Yt(:), [ones(par.N*(par.T-1),1), Yt_1(:)])
            
            % real data
            fprintf('Real data\n');
            Yt = data.Y(:,2:end);
            Yt_1 = data.Y(:,1:end-1);
            regress(Yt(:), [ones(par.N*(par.T-1),1), Yt_1(:)])
        end
        
        function par = update_par_mcmc_from_globals(par)
            global Resqinit_eps Resqinit Resqinit_e0 b1_eps bL_eps b1 bL b1_e0 bL_e0 
            par.eps   = abb_par(Resqinit_eps, b1_eps, bL_eps); 
            par.eta_0 = abb_par(Resqinit_e0, b1_e0, bL_e0); 
            par.eta_t = abb_par(Resqinit, b1, bL); 
        end
        
        function fval = density_wrapper(Matdraw)
            [data,par] = abb2017.repackage_globals(); 
            par = abb2017.update_par_mcmc_from_globals(par); 
            data = abb2017.precompute(data,par); 
            par = abb2017.precompute_par(par); 
            fval = abb2017.density(Matdraw, data, par); 
        end
        
        function par = precompute_par(par)
            % useful for indexing in the hermite basis function computation
            par.precomputed.ii1_K1K2 = kron(1:par.K1, ones(1,par.K2));
            par.precomputed.ii2_K1K2 = kron(ones(1,par.K1), 1:par.K2);
        end
        
        function [data] = precompute(data, par)
            % abb2017.precompute(): age-related precomputations of hermite
            % basis functions. 
            % 
            % Precomputes the following fields: 
            %   data.precomputed_hermite.age_t: MatAGE_t originally
            %   data.precomputed_hermite.age_0: MatAGE1 originally
            

            N = par.N; 
            T = par.T; 
            data.precomputed_hermite = struct(); 
            data.precomputed_hermite.age_t = cell(T,1); 
            
            % 1. full age matrix 
            age = abb2017.standardize_age(data.AGE(:), par); 
            data.precomputed_hermite.age = abb2017.hermite_basis_univariate(age, par.K2); 
            
            % 2. year-by-year
            for t=1:T
                
                % a. find age 
                
                % whether to use code implementing a tiny error from the
                % original code used in the paper (virtually no effect on
                % final results)
                original = false;
                if ~( original && isfield(par, 'Nsim') && par.Nsim>1 )
                    % i. do it right 
                    age = abb2017.standardize_age(data.AGE(:,t), par);
                    
                else
                    % ii. backwards compatible
                    %     the original code for the simulation scrambles
                    %     age in a way that is not compatible with
                    %     recomputation using data.AGE. 
                    if t==1, fprintf('Scrambling hermite basis function\n'); end; 
                    
                    % infer original #obs. (before repmatting) 
                    Npre = par.N/par.Nsim; 
                    assert(Npre == round(Npre), 'Scrambling does not work.'); 
                    
                    % verify hat AGE is copied as this code assumed
                    assert(all(all( repmat(data.AGE(1:Npre, :), par.Nsim,1) == data.AGE )), 'Data is not copied down. Original simulation code cannot be reconstructed, use a global to access MatAGE_t from the abb2017.simulate() function. '); 
                    
                    % construct scrambled in the same way as ABB 
                    AGE = data.AGE(1:Npre,:); 
                    AGE = reshape(repmat(AGE(:), par.Nsim, 1), Npre*par.Nsim, par.T); 
                    age = abb2017.standardize_age(AGE(:,t), par); 
               
                end
                
                % b. precompute 
                data.precomputed_hermite.age_t{t} = abb2017.hermite_basis_univariate(age, par.K4);
            end
            
            % 3. first year 
            age0_std = abb2017.standardize_age(data.AGE(:,1), par); 
            data.precomputed_hermite.age_0 = abb2017.hermite_basis_univariate(age0_std, par.K3); 
        end
        
        function age_std = standardize_age(age, par)
            age_std = (age-par.meanAGE)/par.stdAGE; 
        end
        
        function inc_std = standardize_inc(inc, par)
            inc_std = (inc-par.meanY)/par.stdY; 
        end
                
        function fval = density(Matdraw,data,par)
            [N,T] = size(data.Y); 
            
            
            % ------ eps density ------
            
            Vect = data.Y-Matdraw; 
            Vect = Vect(:); 
            dens = abb2017.inner_density(data.precomputed_hermite.age, Vect, par.eps.Param, par.eps.b1, par.eps.bL, par.Vectau);
            denstot = prod( reshape(dens,N,T), 2);

            
            % ------ eta density -----
            
            dens2=zeros(N,T);
            
            
            % model for initial eta_0 
           
            dens2(:,1) = abb2017.inner_density(data.precomputed_hermite.age_0, Matdraw(:,1), par.eta_0.Param, par.eta_0.b1, par.eta_0.bL, par.Vectau);
            
            
            % iterate for eta_t 
            
            for tt=1:T-1 % NB! think of this as tt=2:T and be mindful of indexing
                
                % --- 1. explanatory variables matrix ---
                
                % 1.a standardize
                Mat_t_std = abb2017.standardize_inc(Matdraw(:,tt), par);   % "Y_lag" since tt and not tt+1
                Herm_inc = abb2017.hermite_basis_univariate(Mat_t_std, par.K1);  
                Herm_age = data.precomputed_hermite.age_t{tt+1}; % to create manually -> %Herm_age = abb2017.hermite_basis_univariate(Age_t_std, K2);
                Mat = abb2017.product_basis(Herm_inc, Herm_age);
                
                % --- 2. density ---
                
                dens2(:,tt+1) = abb2017.inner_density(Mat, Matdraw(:,tt+1), par.eta_t.Param, par.eta_t.b1, par.eta_t.bL, par.Vectau);
                
            end
            
            % take product over T
            dens2tot = prod( reshape(dens2,N,T), 2);
            
            fval=denstot.*dens2tot;
            
        end
        
        function dens = inner_density(Xmat, Yvec, Param_mat, par_b1, par_bL, Vectau)
            % Computes the full likelihood function
            % - Yvec can fall in various intervals demarked by points A_jtau
            % - bounds for interval (Vectau(jtau); Vectau(jtau+1)] is constructed
            %   by A_jtau=Xmat*Param_mat(:,jtau)
            %
            % NOTE: if A_jtau is non-monotonic, then an individual can belong to
            % multiple groups and the likelihood will be added for all groups the
            % individual falls in.
            
            Ntau = numel(Vectau);
            Ny = numel(Yvec);
            assert(isvector(Yvec));
            assert(~isempty(par_bL), 'input, par_bL, was empty!'); 
            
            dens = zeros(Ny,1);
            
            % 1. interior
            
            for jtau=1:Ntau-1
                
                % 1.a bounds
                A1 = Xmat*Param_mat(:,jtau);
                A2 = Xmat*Param_mat(:,jtau+1);
                I = (Yvec>A1) & (Yvec<=A2);
                
                % 1.b density
                diff_tau = Vectau(jtau+1)-Vectau(jtau);
                diff_A = (A2 - A1);
                dens(I) = dens(I) + diff_tau./diff_A(I);
                
            end
            
            % 2. extreme intervals
            
            % 2.a low
            A_low  = Xmat*Param_mat(:,1);
            I = Yvec<=A_low;
            dens(I) = dens(I) + Vectau(1)*par_b1*exp(par_b1*(Yvec(I)-A_low(I)));
            
            % 2.b high
            A_high = Xmat*Param_mat(:,Ntau);
            I = Yvec>(A_high);
            dens(I) = dens(I) + (1-Vectau(Ntau))*par_bL*exp(-par_bL*(Yvec(I)-A_high(I)));
            
        end
        
        function [Y, eta, eps] = simulate(data, par, pmin_vec, pmax_vec) 
            
            assert(isfield(data, 'V_draw_eta'), 'Missing field "V_draw_eta" from data'); 
            assert(isfield(data, 'V_draw_eps'), 'Missing field "V_draw_eps" from data'); 
            
            N = par.N; 
            T = par.T; 
            assert( all(size(data.V_draw_eta) == [N,T]), 'dimension failure'); 
            
            % ------------------- draw eta (persistent component) ---------------------
            
            eta=zeros(N,T);
            
            % eta_0 (initial conditions)
            V_draw = data.V_draw_eta(:,1);
            eta(:,1) = abb2017.draw_eta_t(V_draw, data.precomputed_hermite.age_0, par.eta_0.Param, par.Vectau, par.eta_0.b1, par.eta_0.bL);
            
            % eta_t (recursive version)
            for tt=1:T-1 %
                
                % 1. set up
                var1 = (eta(:,tt)-par.meanY)/par.stdY;
                var2 = (data.AGE(:,tt+1)-par.meanAGE)/par.stdAGE;
                Mat = abb2017.hermite_basis_bivariate(var1,var2,par.K1,par.K2);
                
                % Alternative, faster
                %H_y = abb2017.hermite_basis_univariate(var1, par.K1);
                %Mat2 = abb2017.product_basis(H_y, data.precomputed_hermite.age_t{tt+1}); 
                %assert(mean(mean( abs(Mat - Mat2) ))<1e-8) ;
                
                %V_draw=unifrnd(0,1,N,1);
                V_draw = data.V_draw_eta(:,tt+1); 
                
                % censoring thresholds
                pmin = pmin_vec(tt+1);
                pmax = pmax_vec(tt+1); 
                eta(:,tt+1) = abb2017.draw_eta_t(V_draw, Mat, par.eta_t.Param, par.Vectau, par.eta_t.b1, par.eta_t.bL, pmin, pmax);
                
            end
            
            
            
            % ------------------- draw eps (transitory component) ---------------------
            
            eps=zeros(N,T);
            
            for tt=1:T
                % Proposal, eta_0
                V_draw = data.V_draw_eps(:,tt); 
                
                MatAGE_t_ptr = data.precomputed_hermite.age_t{tt};
                
                eps(:,tt) = abb2017.draw_eta_t(V_draw, MatAGE_t_ptr, par.eps.Param, par.Vectau, par.eps.b1, par.eps.bL);
                
            end
            
            Y=eta+eps;
        end
        
        function eta_draws = draw_eta_t(V_draw, Xmat, Params, Vectau, par_b1, par_bL, cut_min, cut_max)
            % cut_{min,max}: censoring thresholds: simulated eta will be 
            % cut at that value. set to -inf, +inf if you do not want. 
            % ABB (2017) use 3*{min/max}(Y) 
            
            N = numel(V_draw);
            Ntau = numel(Vectau);
            
            assert(isvector(V_draw));
            assert(~isempty(Xmat), 'No Xmat input!'); 
            assert(all(size(Params) == [size(Xmat,2),Ntau]));
            assert(size(Xmat,1) == size(V_draw,1));
            
            eta_draws = zeros(N,1);
            
            % first quantile
            I = V_draw<=Vectau(1);
            eta_draws(I) = Xmat(I,:)*Params(:,1);
            
            % middle
            for jtau=2:Ntau
                I = (V_draw>Vectau(jtau-1)) & (V_draw<=Vectau(jtau));
                A1 = Xmat*Params(:,jtau);
                A0 = Xmat*Params(:,jtau-1);
                dA_dtau = (A1(I)-A0(I))/(Vectau(jtau)-Vectau(jtau-1));
                eta_draws(I,1)=eta_draws(I,1) ...
                    +( dA_dtau.*(V_draw(I)-Vectau(jtau-1)) + A0(I) );
            end
            
            %Last quantile.
            I = (V_draw>Vectau(Ntau));
            eta_draws(I,1)=eta_draws(I,1)+(Xmat(I,:)*Params(:,Ntau));
            
            % Laplace tails
            Ihi = (V_draw>Vectau(Ntau));
            eta_draws(Ihi,1)=eta_draws(Ihi,1) -(1/par_bL *log((1-V_draw(Ihi))/(1-Vectau(Ntau))));
            
            Ilo = (V_draw<=Vectau(1));
            eta_draws(Ilo,1)=eta_draws(Ilo,1) +(1/par_b1 *log(V_draw(Ilo)/Vectau(1)));
            
            % censor at top/bottom
            if nargin>=7 
                if ~isnan(cut_min)
                    eta_draws = max(eta_draws, cut_min); 
                end
                if ~isnan(cut_max) 
                    eta_draws = min(eta_draws, cut_max); 
                end
                %eta_draws = min(max(eta_draws, cut_min), cut_max);
            end
            
        end
        
        function [b1,bL] = update_b_par(Yvec, Xmat, Param, Ntau)
            % updating of the parameters governing the density behavior in
            % the tails (first and last sub-intervals of [0;1]) is
            % determined by simple moment equations directly. 
            
            assert(size(Yvec,2)==1); 
            assert(size(Param,2) == Ntau); 
            assert(size(Xmat, 1) == size(Yvec, 1)); 
            
            Vect1 = Yvec - Xmat*Param(:,1); 
            I = Vect1<=0; 
            b1 = -sum(I) / (Vect1'*I);
            
            Vect2 = Yvec - Xmat*Param(:,Ntau); 
            I = Vect2>=0; 
            bL = sum(I)/(Vect2'*I);
            
        end
        
        function H = product_basis(H1, H2, par)
            [N ,K1] = size(H1); 
            [N2,K2] = size(H2); 
            assert(N==N2); 
            
            H = nan(N, K1*K2); 
            
            % Fast indexing 
            if nargin>=3 && isfield(par, 'precomputed') && isfield(par.precomputed, 'ii1_K1K2') && isfield(par.precomputed, 'ii2_K1K2')
                % if precomputed
                ii1 = par.precomputed.ii1_K1K2;
                ii2 = par.precomputed.ii2_K1K2;
            else
                % do it here 
                ii1 = kron(1:K1, ones(1,K2)); 
                ii2 = kron(ones(1,K1), 1:K2); 
            end
                        
            H = H1(:,ii1) .* H2(:,ii2); 
            
            % slow but readable
            % i = 1; 
            % for k1=1:K1
            %     for k2=1:K2
            %         H(:,i) = H1(:,k1).*H2(:,k2); 
            %         i = i+1; 
            %     end
            % end
            
        end
        
        function H = hermite_basis_bivariate(x1,x2,K1,K2)
            % 1. precompute both 
            H1 = abb2017.hermite_basis_univariate(x1, K1); 
            H2 = abb2017.hermite_basis_univariate(x2, K2); 
            
            % 2. multiply 
            H = abb2017.product_basis(H1,H2); 
        end
        
        function H = hermite_basis_univariate(x, K)
            % x = vector 
            % K = order of polynomium - 1 (because 0 also is included
            N = size(x,1); 
            assert(size(x,2) == 1, 'x not univariate');
            
            H = nan(N,K+1); 
            for k=0:K
                i = k+1; 
                H(:,i) = hermite(k, x); 
            end
            
        end
        
        
        function [Nu_vec2, likelihood_vec2, acceptrate] = update_chain(data, par, Nu_vec1, likelihood_vec1)
            % initialize draws: same as last period
            assert(all(size(Nu_vec1) == [par.N, par.T]), 'Nu_chain dimensions failure') ;
            assert(all(size(likelihood_vec1) == [par.N, 1]), 'Obj_chain dimensions failure') ;
            
            Matdraw = Nu_vec1;
            Nu_vec2 = nan(size(Nu_vec1)); 
            acceptrate = nan(par.T,1); 
            
            for t=1:par.T
                
                if t==1
                    % intialization: from previous j
                    oldObj = likelihood_vec1;
                else
                    % iterate on current value
                    oldObj = likelihood_vec2;
                end
                
                % new draw
                Matdraw(:,t)=Nu_vec1(:,t)+sqrt(par.var_prop(t))*randn(par.N,1);
                
                % compute new density
                %newObj=postr_QRMCMC_age_hermite(Matdraw);
                %newObj = abb2017.density_wrapper(Matdraw);
                newObj = abb2017.density(Matdraw, data, par);
                
                
                % accept/reject draws
                r = min(newObj./oldObj, 1); % probabilities must be <= 1
                Iaccept=rand(par.N,1)<=r;
                
                % obj
                likelihood_vec2 = oldObj; 
                likelihood_vec2(Iaccept)  = newObj(Iaccept);
                
                % nu
                Nu_vec2(:,t) = Nu_vec1(:,t); 
                Nu_vec2(Iaccept,t)  = Matdraw(Iaccept,t);
                
                % save
                Matdraw(:,t)=Nu_vec2(:,t);
                acceptrate(t) = mean(Iaccept);
                
            end
        end
        
        
        function [par, Matdraw] = m_step(Nu_draws,data,par)
            
            N = par.N; 
            T = par.T; 
            
            % Update explanatory variables 
            
            %Last draws of the chain will be the fixed associated with our data.
            Matdraw = Nu_draws;
            
            Matdraw_t     = reshape(Matdraw(:,2:T),   N*(T-1), 1); % y-var for eta_t
            Matdraw_tot   = reshape(Matdraw,          N*T,     1); % part of y-var for eps
            
            % form cartesian basis functions in {lagged draw}*{age}
            %Matdraw_lag=nan(N*(T-1), (K1+1)*(K2+1));
            Matdraw_t_lag = reshape(Matdraw(:,1:T-1), N*(T-1), 1);
            var1 = (Matdraw_t_lag-par.meanY)/par.stdY;
            age  = reshape(data.AGE(:,2:end), N*(T-1), 1); % equivalent to MatAGE_t_tot when Mdraws==1 
            var2 = (age-par.meanAGE)/par.stdAGE;
            Matdraw_lag = abb2017.hermite_basis_bivariate(var1,var2,par.K1,par.K2); % x-var for eta_t
            
            % these really do assume Mdraws == 1 
            MatAGE_tot = data.precomputed_hermite.age; 
            MatAGE1_tot = data.precomputed_hermite.age_0; 
            Ytot_t = data.Y(:); 
            
            % ----------------------------- M step --------------------------------
            
            Ntau = numel(par.Vectau); 
            
            for jtau=1:Ntau
                
                tau=par.Vectau(jtau);
                
                % 1. eta_0: y=first inc, x = first age
                par.eta_0.Param(:,jtau) = rq(MatAGE1_tot, Matdraw(:,1), tau);
                
                % 2. eta_t: y=inc, x=lag_inc
                par.eta_t.Param(:,jtau) = rq(Matdraw_lag, Matdraw_t, tau);
                
                % 3. eps:   y=inc-draws x=age
                par.eps.Param(:,jtau) = rq(MatAGE_tot, Ytot_t-Matdraw_tot, tau);
                
            end

            % Normalization
            par.eps.Param=par.eps.Param-mean(par.eps.Param')'*ones(1,Ntau);
            par.eps.Param(1,:)=par.eps.Param(1,:)-((1-par.Vectau(Ntau))/par.eps.bL-par.Vectau(1)/par.eps.b1)*ones(1,Ntau);
            
            
            % update Laplace parameters
            [par.eps.b1,par.eps.bL]     = abb2017.update_b_par(Ytot_t-Matdraw_tot, MatAGE_tot, par.eps.Param, Ntau); 
            [par.eta_0.b1,par.eta_0.bL] = abb2017.update_b_par(Matdraw(:,1), MatAGE1_tot,      par.eta_0.Param, Ntau); 
            [par.eta_t.b1,par.eta_t.bL]  = abb2017.update_b_par(Matdraw_t, Matdraw_lag,         par.eta_t.Param, Ntau); 
            
        end
        
        
        
        function quick_persistence_computation_mat(Nu_draws, data, par)
            
            N = par.N; 
            T = par.T; 
            K1 = par.K1; 
            K2 = par.K2; 
            Ntau = numel(par.Vectau); 
                
            Vect=Nu_draws(:,1:T-1);          % lagged Y
            Vect=quantile(Vect(:),par.Vectau);  % take Ntau quantiles of lagged Y
            age_ref=par.meanAGE;                % age at which to evaluate
        
            % matrix of regressors
            Mat=zeros(Ntau,K2+1);
            for kk1=1:K1
                for kk2=0:K2
                    Mat=[Mat kk1*hermite(kk1-1,(Vect(:)-par.meanY)/par.stdY)./par.stdY.*hermite(kk2,(age_ref-par.meanAGE)/par.stdAGE)];
                end
            end
            
            fprintf('rho(perc. of Y_t-1, perc. of Y_t):\n');
            persistence = Mat*par.eta_t.Param; 
            disp(persistence); 
        end
        
        
        
        function quick_persistence_computation_avg(Nu_draws, data, par)
            
            N = par.N; 
            T = par.T; 
            K1 = par.K1; 
            K2 = par.K2; 
            Ntau = numel(par.Vectau); 
            
            % quick computation of persistence
            Vect=Nu_draws(:,1:T-1);
        
            Mat=zeros(N*(T-1),K2+1);
            AGE_t = reshape(data.AGE(:,2:end), N*(T-1), 1); 
            for kk1=1:K1
                for kk2=0:K2
                    Mat=[Mat kk1*hermite(kk1-1,(Vect(:)-par.meanY)/par.stdY)./par.stdY.*hermite(kk2,(AGE_t-par.meanAGE)/par.stdAGE)];
                end
            end
        
            meanVec = mean(Mat*par.eta_t.Param); 
            
            fprintf('%14s', 'tau = '); 
            for jtau=1:Ntau 
                fprintf('%8.4f ', par.Vectau(jtau))
            end
            fprintf('\n%14s', 'Mean rho = '); 
            for jtau=1:numel(par.Vectau) 
                fprintf('%8.4f ', meanVec(jtau)); 
            end
            fprintf('\n'); 
            
        end
        
        function estimates_chain = update_estimates_chain(par, iter)
            estimates_chain = par.estimates_chain; 
            estimates_chain.eta_t(iter) = par.eta_t; 
            estimates_chain.eta_0(iter) = par.eta_0; 
            estimates_chain.eps(iter)   = par.eps; 
        end
        
        function mean_param = mean_param_from_chain(estimates_chain_par, parName, ii)
            % EXAMPLE: 
            % abb2017.mean_param_from_chain(par.estimates_chain.eta_t,'Param', 1:draws)
            assert(isprop(estimates_chain_par(1), parName), 'Parameter "%s" not found in estimates_chain(1).', parName);
            maxiter = numel(estimates_chain_par); 
            if nargin<3
                ii = round(ii/2):maxiter; 
            end
            
            % look for emtpy 
            isEmpty = cellfun(@isempty, {estimates_chain_par(ii).Param}); 
            assert(~any(isEmpty), 'Empty parameters found in estimates_chain for requested range, ii. Maybe you are not done estimating? '); 
            
            siz = size(estimates_chain_par(1).(parName)); 
            pars = [estimates_chain_par(ii).(parName)];  % pipe all values out, stack vertically
            pars = reshape(pars, siz(1), siz(2), numel(ii)); 
            mean_param = mean(pars, 3);

        end
        
        function singlepar = mean_all_param_from_chain(chain, ii)
            % mean_all_param_from_chain: compute average parameters from a
            % Monte Carlo chain of estimates
            %
            % INPUTS: 
            %   chain: a list of parameter structs 
            %   ii: list of integers; the estimates to use. E.g. 50:100 to
            %   drop the first 50. 
            props = {'Param','b1','bL'}; 
            args = cell(size(props)); 
            for i=1:numel(props)
                var = props{i}; 
                assert(all(isprop(chain, var)), 'Property, "%s", not found in chain.', var);
                args{i} = abb2017.mean_param_from_chain(chain, var, ii); 
            end
            singlepar = abb_par(args{:}); 
        end
        
        function par = assign_par_to_mean_of_chain(par, ii)
            par.eta_t = abb2017.mean_all_param_from_chain(par.estimates_chain.eta_t, ii); 
            par.eta_0 = abb2017.mean_all_param_from_chain(par.estimates_chain.eta_0, ii); 
            par.eps   = abb2017.mean_all_param_from_chain(par.estimates_chain.eps, ii); 
        end

        function update_globals_from_par_mcmc(par)
            % legacy function no longer in use. However, it shows how to
            % translate from the original code variable names to the par
            % struct counterparts. 
            global Resqinit_eps Resqinit Resqinit_e0 b1_eps bL_eps b1 bL b1_e0 bL_e0 
            Resqinit_eps = par.eps.Param; 
            Resqinit     = par.eta_t.Param; 
            Resqinit_e0  = par.eta_0.Param; 
            b1_eps       = par.eps.b1; 
            bL_eps       = par.eps.bL; 
            b1_e0        = par.eta_0.b1; 
            bL_e0        = par.eta_0.bL; 
            b1           = par.eta_t.b1; 
            bL           = par.eta_t.bL; 
        end

        function plot_moment_comparisons(Ysim, Ydat, first_age, DOSAVE)
            
            if isempty(first_age)
                first_age = 30; 
            end
            
            [N,T] = size(Ysim); 
            assert(all( size(Ydat) == [N,T] ), 'Simulated and data matrices do not conform.');
            
            aa = (1:T) + first_age-1; 
            
            % Compare moments
            
            figure(1);
            f = @(x) mean(x)';
            plot(aa,f(Ysim),'-o', aa,f(Ydat),'-x'); legend('Simulation','Data'); ylabel('Mean');
            if DOSAVE
                saveas(gcf, '../output/abb_sim_vs_dat_mean.pdf'); 
            end
            
            figure(2);
            f = @(x) var(x)';
            plot(aa,f(Ysim),'-o', aa,f(Ydat),'-x'); legend('Simulation','Data'); ylabel('Variance');
            if DOSAVE
                saveas(gcf, '../output/abb_sim_vs_dat_var.pdf');
            end
            
            figure(3);
            f = @(x) skewness(x)';
            plot(aa,f(Ysim),'-o', aa,f(Ydat),'-x'); legend('Simulation','Data'); ylabel('Skewness');
            if DOSAVE
                saveas(gcf, '../output/abb_sim_vs_dat_skew.pdf');
            end
            
            figure(4);
            f = @(x) kurtosis(x)';
            plot(aa,f(Ysim),'-o', aa,f(Ydat),'-x'); legend('Simulation','Data'); ylabel('Kurtosis');
            if DOSAVE 
                saveas(gcf, '../output/abb_sim_vs_dat_kurt.pdf');
            end
        end
        
        function print_avg_acceptrate(acceptrate)
            T = size(acceptrate, 1); 
            fprintf('Accept rates: ');
            for t=1:T
                fprintf('t=%d: %5.2f%%  ', t, 100.0*mean(acceptrate(t,:)));
            end
            fprintf('\n');
        end
        
        function print_parameters(par)
            fprintf('--- eta_0 ---\n'); 
            par.eta_0.print();
            fprintf('--- eta_t ---\n'); 
            par.eta_t.print();
            fprintf('--- eps ---\n'); 
            par.eps.print();
        end
        
    end
end