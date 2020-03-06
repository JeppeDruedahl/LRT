%% ####################### ABB 2017: Simulate #############################
% 
% Code for Arellano, Blundell and Bonhomme (2017; ECTA), 
%   "Earnings and Consumption Dynamics: A Nonlinear Panel Data Framework"
%
% This code is based on the original code accompanying the paper but
% modified by Anders Munk-Nielsen to work for general income datasets. 
%
% CORE DATA: 
%   Y: N*T matrix of log income observations. Read in from the csv file, 
%      ../data/logY_p1.csv: A CSV file with N rows (individuals) and T
%      columns (ages). 
%   AGE: N*T matrix of household ages. 
% 

close all; 
clear all
clc;

%% ############################## SETTINGS ################################

% Number of ages to include in estimation sample (randomly draw start age)
T = 15;

% grid of taus over [0;1] 
Ntau = 11;
Vectau = (1/(Ntau+1):1/(Ntau+1):Ntau/(Ntau+1))'; 

% Hermite polynomial degrees
K1 = 3; % eta_t 
K2 = 2; % eta_t 
K3 = 2; % eta_0 
K4 = 2; % eps 
KPolyAgeResid=5; % residualization of income against age polynomial

% variance Random Walk proposals
var_prop = [0.08; repmat(0.03, T-2, 1); 0.05]; 

% Complexity 
maxiter = 500; % Markov chain iterations (default=500) 
draws = 200;   % Number of draws within the chain

% Technical 
rng('default')
USEWAITBAR = true;  % show progress with a nice graphical waiting bar

%% ############################### READ ###################################

% 1. log income data file 
Y = csvread('../data/logY_p100.csv');
%Y = Y(1:100, :); 

% 2. age matrix (assuming the N*T Y-matrix has ages 30, 31, ..., 59. 
tt = (1:size(Y,2)) + 29; % assumes first age is 30; last is 30+size(Y,2)-1
AGE = repmat(tt, size(Y,1), 1);

% 3. put into structs 
data = abb2017.initialize(Y, AGE, T, K1, K2, K3, K4, KPolyAgeResid, Vectau, var_prop, draws, maxiter); 


%% ############################# SIMULATE #################################

% 1. simulate data with same N 
Nsim = size(Y, 1); 
ages_sim = 30:59; 
par = read_par_estimates(); 
par.T = 30; 
Ysim = abb2017.simulate_data(Nsim, ages_sim, data, par); 

% 2. write to disk 
csvwrite('../data/abb_sim.csv', Ysim); 
fprintf('Simulation written to disk.\n'); 


%% ############################### PLOT ###################################

% 1. do a small sanity check 
data = abb2017.initialize(Y, AGE, T, K1, K2, K3, K4, KPolyAgeResid, Vectau, var_prop, draws, maxiter); 
[par.N, par.T] = size(data.Y); 
abb2017.simulate_and_compare_AR1_regression(data, par); 

% 2. compare moments on full sample 
first_age = 30; % t=1 corresponds to age 30 
DOSAVE = false; 
abb2017.plot_moment_comparisons(Ysim, Y, first_age, DOSAVE)

