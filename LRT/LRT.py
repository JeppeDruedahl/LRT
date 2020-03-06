import time
import copy
import pickle 
import itertools as it

import numpy as np
import pandas as pd

from sklearn import tree

class ParStruct:

    def __init__(self):

        # a. general
        self.T = 30 # time periods
        self.N = 100000 # individuals
        
        self.agemin = 30         
        self.agemax = 59
        
        self.simN = 200000
        
        # b. moments
        self.mainmoments_varname_str_list = ['dlogY']
    
        # recent earnings
        self.RE_lags = 5
        self.RE_min_obs = 3
        
        # age groups
        self.age_grps = dict()        
        self.age_grps['all'] = []
        for i in range(self.T): 
            self.age_grps['all'].append(i+self.agemin)
        self.age_grps['detail'] = [[30,34],[35,39],[40,44],
                                   [45,49],[50,54],[55,59]]
                                   
        # percentiles groups
        self.perc_grps = dict()

        self.perc_grps['all'] = []
        for i, j in zip(np.linspace(0,99,100),np.linspace(1,100,100)):
            self.perc_grps['all'].append([i,j])
        
        self.perc_grps['RE'] = [[0,1],[1,10],[10,20],[20,30],
                                [30,40],[40,50],[50,60],[60,70],
                                [70,80],[80,90],[90,95],[95,99],[99,100]]                 
        
        # leads and lags
        self.K_leads = [1,5,10]   
        self.K_dlogRE_leads = [5]             
        self.K_autocorr = [1,2,5,10]        
        
        # cov(y,dy): which s,t,k to be used 
        self.cov_YdY_tt = [10,15,20]
        self.cov_YdY_ss = [5,10,15]
        self.cov_YdY_kk = np.arange(1,16) # to but not including end-point

        # midpoints of percentiles
        self.perc_x = dict()
        for perc_grp_str in ['all','RE']:
            perc_grp = self.perc_grps[perc_grp_str]
            self.perc_x[perc_grp_str] = np.zeros(len(perc_grp))
            x = self.perc_x[perc_grp_str]
            for ipercs, perc_bounds in enumerate(perc_grp):
                x[ipercs] = np.mean(np.array(perc_bounds))
                
class modelStruct:

    def __init__(self):
        
        self.tree = []
        self.uniqueG = []
        self.ypred_G = []
        self.trans = []
        self.trans_obs = []        
        self.G = []
        self.I = []

class dataStruct: 
    def __init__(self):
        pass

def mse_scorer(estimator,x,y): 

    ypred = estimator.predict(x)
    mse = np.mean((ypred-y)**2)
    return mse 

def get_features_and_targets(par,data,t): 

    [T,N] = data.logY.shape 

    # a. dimensions and allocations
    if hasattr(par, 'use_exogenous_features'): 
        USE_EXO_FEATURES = par.use_exogenous_features
    else: 
        USE_EXO_FEATURES = False

    # i. allocate x
    tlags = np.arange(np.fmax(t-par.k,0),t+1)
    if USE_EXO_FEATURES: 
        sh = data.exogenous_features[t].shape
        if len(sh) == 1: 
            dim_exo = 1
        else: 
            dim_exo = sh[1]
        x = np.empty((N,tlags.size + dim_exo)) 
    else: 
        x = np.empty((N,tlags.size))   

    # ii. allocate y 
    tleads = np.arange(t,np.fmin(t+par.k_lead,T-1)+1)        
    y = np.empty((N,tleads.size))

    # b. create x

    # i. find lags 
    for i, tlag in enumerate(tlags):
        x[:,i] = data.logY[tlag,:]

    # ii. fetch exogenous features 
    if USE_EXO_FEATURES: 
        if dim_exo == 1:
            x[:,tlags.size] = data.exogenous_features[t]
        else: 
            x[:,tlags.size:] = data.exogenous_features[t]

    # c. create y
    for i, tlead in enumerate(tleads):
        y[:,i] = data.logY[tlead,:]

    return x,y

def estimate(par,data,name,color='black'):
    
    # 1. initialize result struct
    model = modelStruct()
    model.type = 'ML'
    model.name = name
    model.par = copy.deepcopy(par)
    model.color = color
       
    Tdata, _Ndata = data.logY.shape
    model.num_leafs = np.zeros(Tdata,dtype=np.int)

    # 2. build trees and construct transition matrices    
    for t in range(par.T):
        
        printfac = 5
        if t%printfac == 0:
            if t == 0:
                print("t = {}".format(t))
            else:
                print("t = {:2d} ({:3.2f} secs per period)"
                    .format(t,(time.time()-t1)/printfac))
            t1 = time.time()

        # a. data
        t0 = time.time()
        
        # i. find lags/leads 
        x,y = get_features_and_targets(par, data, t)
            
        # ii. handle missings
        ImissX = np.amax(np.isnan(x),axis=1)
        ImissY = np.amax(np.isnan(y),axis=1)
        Ikeep = (ImissX == 0) & (ImissY == 0)

        x = x[Ikeep]
        y = y[Ikeep]

        if printfac == 1:
            print(" data:   {:3.2f} secs".format(time.time()-t0))
            
        # b. construct tree
        t0 = time.time()
        
        treenow = tree.DecisionTreeRegressor(max_depth=par.depth,min_samples_leaf=100)
        treenow.fit(x,y)
        model.tree.append(treenow)

        if printfac == 1:
            print(" tree:   {:3.2f} secs".format(time.time()-t0))
            
        # c. find groups
        
        t0 = time.time()
        
        # i. find
        G = -1*np.ones((len(Ikeep),), dtype=int)
        G[Ikeep] = treenow.apply(x)
        model.G.append(G.astype(np.int))
        
        # ii. unique
        uniqueG = np.unique(G)        
        J = uniqueG >= 0 # remove missings
        uniqueG = uniqueG[J]
        model.uniqueG.append(uniqueG.astype(np.int))

        # iii. number of leafs        
        num_leafs = len(uniqueG)
        model.num_leafs[t] = num_leafs
        
        if printfac == 1:
            print(" groups: {:3.2f} secs".format(time.time()-t0))
            
        # d. find prediction for each group
        t0 = time.time()
        
        ypred_all = treenow.predict(x)
        
        # i. only the first elements
        _ydmim1, ydim2 = y.shape
        if ydim2 == 1:
            ypred = ypred_all
        else:
            ypred = ypred_all[:,0] 
        
        # ii. find unique prediction within groups
        ypred_G = np.empty(len(uniqueG))                
        for j, Gnow in enumerate(uniqueG):
            J = G[Ikeep] == Gnow
            ypred_G[j] = ypred[J][0] # ypred[J] is unique
        model.ypred_G.append(ypred_G)
        
        if printfac == 1:
            print(" pred:   {:3.2f} secs".format(time.time()-t0))
            
        # e. transition matrix
        t0 = time.time()
        
        if t == 0:            
            model.trans.append([])
            model.trans_obs.append([])            
            continue        
        
        trans = np.empty((model.num_leafs[t-1],model.num_leafs[t]))   
        trans_obs = np.empty((model.num_leafs[t-1],model.num_leafs[t]))           
        maxG = np.amax(model.uniqueG[t])
        for i in range(model.num_leafs[t-1]): 
            
            Q = (model.G[t-1] == model.uniqueG[t-1][i]) & (model.G[t] > -1)

            # i. count transitions from G_t-1==i to all G_t 
            countG = np.bincount(model.G[t][Q], minlength=maxG+1)            
            trans[i,:] = countG[model.uniqueG[t]]
            
            # ii. normalize 
            Nlag = np.sum(trans[i,:])
            if not Nlag > 0: 

                trans_obs[i,:] = 0
                trans[i,:] = 0.0    

            else:
                
                trans_obs[i,:] = np.copy(trans[i,:])
                trans[i,:] = trans[i,:] / Nlag 

        model.trans.append(trans)
        model.trans_obs.append(trans_obs)

        if printfac == 1:
            print(" trans:  {:3.2f} secs".format(time.time()-t0))
                    
    return model        

def simulate(par,model,data_in,t0=0,tend='par.T',seed=1917,rng_state=None):

    T,N = data_in.logY.shape
    assert T == par.T, 'T-dimension disagreement between data_in and par'
    assert N == par.N, 'N-dimension disagreement between data_in and par'

    if not seed == None: 
        np.random.seed(seed)
    else:
        np.random.set_state(rng_state)

    # a. default end 
    if tend=='par.T': 
        tend = par.T

    # b. allocate        
    data = dataStruct()
            
    data.logY = np.nan*np.ones((par.T,par.simN))
    data.G = np.zeros((par.T,par.simN),dtype=np.int)
    data.age = np.zeros((par.T,par.simN)) 

    aa = np.arange(par.agemin,par.agemax+1)  
    data.age = data.age + np.tile(aa[:,np.newaxis],[1,par.simN])

    # additional exogenous features
    if hasattr(par, 'use_exogenous_features'): 
        USE_EXO_FEATURES = par.use_exogenous_features
    else:
        USE_EXO_FEATURES = False

    # c. loop
    first = 1
    for t in range(t0,tend):
                
        # i. first period
        if first == 1:            
            
            first = 0
            
            # o. initial groups too draw from
            tlags = np.arange(np.fmax(t-par.k,0),t+1) 
            if USE_EXO_FEATURES: 
                exo = data_in.exogenous_features[t]
                assert exo.shape[0] == par.N , f'exogenous features has {N_exo} elements; par.N is {par.N}'
                if len(exo.shape) == 1: 
                    dim_exo = 1
                    N_exo = exo.shape
                else: 
                    N_exo,dim_exo = exo.shape
            else: 
                dim_exo = 0

            x = np.empty((par.N,tlags.size + dim_exo))        
            for i, tlag in enumerate(tlags):
                x[:,i] = data_in.logY[tlag,:]

            if USE_EXO_FEATURES: 
                if dim_exo == 1:
                    x[:,tlags.size] = exo
                else: 
                    x[:,tlags.size:] = exo

            ImissX = np.amax(np.isnan(x),axis=1)    
            GnonmissX = model.tree[t].apply(x[~ImissX,:].astype(np.float32))    
            
            # oo. initial group 
            if (par.simN == par.N) & (np.any(ImissX) == False): 
                # just use the data as initial period 
                G = GnonmissX
            else: 
                # random group (with replacement); probability proportional to data  
                G = np.random.choice(GnonmissX,par.simN)

            # ooo. convert to 0,1,2,3,...
            for i in range(model.num_leafs[t]):
                I = G == model.uniqueG[t][i]
                data.G[t,I] = i
            
        # ii. all other periods
        else:
            
            # o. groups in current period
            num_leafs = model.num_leafs[t] 
            
            # oo. loop over groups in last period
            for i in range(model.num_leafs[t-1]):
                
                I = data.G[t-1,:] == i
                num = np.sum(I)
                                               
                if num > 0:                                                      
                    pvec = model.trans[t][i,:]
                    data.G[t,I] = np.random.choice(num_leafs,num,p=pvec)   
            
        # iii. earnings
        data.logY[t,:] = model.ypred_G[t][data.G[t,:]]
    
    return data

def compute_MISE_analytically(model,logY,t0=0,t1='T-1',modeltype='LRT'): 
    
    T,N = logY.shape
    if t1=='T-1': 
        t1 = T-1
            
    # a. classify income to discrete gruops 
    
    # LRT 
    if modeltype == 'LRT':
        
        # i. assign input features 
        tlags = np.arange(np.fmax(t0-model.par.k,0),t0+1)  
        x = np.empty((N,tlags.size))        
        for i, tlag in enumerate(tlags):
            x[:,i] = logY[tlag,:]

        G = model.tree[t0].apply(x.astype(np.float32))  
        
        
        # i. convert to base 0
        G_base0 = -1*np.ones(G.shape, dtype=np.int)
        for i,g in enumerate(model.uniqueG[t0]): 
            I = G==g
            if np.any(I): 
                G_base0[I] = i
            
    # PCTBIN 
    elif modeltype == 'PCTBIN':
        y0 = logY[t0,:]
        G_base0 = pd.cut(y0, model.bins[t0], labels=False).astype(int)

    # b. transition matrix from groups in t0 to groups in t1 
    # i. initialize
    prob = model.trans[t0+1]

    # ii. iterate
    for t in range(t0+1,t1): 
        T1 = model.trans[t+1] 
        prob = np.dot(prob, T1) 

    # iii.
    trans_t0_t1 = prob[G_base0, :]

    # c. form error 
    
    # i. tall matrix of errors 
    incs = model.ypred_G[t1]
    err = np.subtract.outer(logY[t1,:], incs)
    sq_err = err**2
    
    # ii. integrate 
    integrated_sq_err = np.sum(sq_err*trans_t0_t1, axis=1)
    
    # iii. avg. across individuals 
    MISE = np.mean(integrated_sq_err)
    
    return MISE
