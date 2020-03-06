import copy
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import scipy 
import copy 

from . import LRT
from . import figfuns

def estimate(par,data,num_bins):

    T,_N = data.logY.shape
    
    # a. setup
    model = LRT.modelStruct() 
    model.name = f'PCTBIN'
    model.name_short = f'PCTBIN_{num_bins}'
    model.par = copy.deepcopy(par)    
    model.num_bins = num_bins
    model.bins = []
    
    # b. construct age-by-age
    for t in range(T): 

        # i. group
        G_t, bins_t = pd.qcut(data.logY[t,:],q=num_bins,labels=False,retbins=True)
        bins_t[0] = -np.inf
        bins_t[-1] = np.inf
        model.G.append(G_t)
        model.bins.append(bins_t) 
        
        # ii. avg. income within bins 
        ypred_G = pd.DataFrame({'y':data.logY[t,:], 'bin':G_t}).groupby('bin').y.mean().values
        model.ypred_G.append(ypred_G)          
        
        # iii. transitions 
        if t == 0: 

            model.trans.append([])

        else: 

            trans = np.zeros((num_bins,num_bins))
            trans_obs = np.zeros((num_bins, num_bins))

            for i_past in range(num_bins): 
                I = model.G[t-1] == i_past
                Nlag = np.sum(I)
                if Nlag > 0: 
                    trans_obs[i_past, :] = np.bincount(model.G[t][I], minlength=num_bins)
                    trans[i_past, :] = trans_obs[i_past, :]/Nlag
                
            model.trans.append(trans)
            model.trans_obs.append(trans_obs)
            
    return model 

def classify(y, bins):

    G = pd.cut(y,bins,labels=False).astype(int)
    return G

def simulate(model,data_in,seed=1917,rng_state=None): 
    
    if not seed == None: 
        np.random.seed(seed)
    else:
        np.random.set_state(rng_state)

    # a. seutp
    T,N = data_in.logY.shape
    num_bins = model.num_bins 

    data = LRT.dataStruct() 
    data.G = np.empty((T,N), dtype=np.int)
    data.logY = np.empty((T,N))
    
    # b. simulate
    for t in range(T): 
        if t == 0: 
            
            data.G[0,:] = classify(data_in.logY[0,:], model.bins[0])
        
        else: 

            for i_past in range(num_bins): 

                I = data.G[t-1,:] == i_past
                num = np.sum(I)
                
                if num > 0:

                    pvec = model.trans[t][i_past,:]
                    data.G[t,I] = np.random.choice(num_bins, num, p=pvec)
        
        data.logY[t,:] = model.ypred_G[t][data.G[t,:]]
    
    return data 