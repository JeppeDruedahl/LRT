import copy
import numpy as np

from . import LRT
from . import moments

def estimate(par,data): 

    # a. setup
    model = LRT.modelStruct() 
    model.name = f'PT'
    model.name_short = f'PT'
    model.par = copy.deepcopy(par)
    data = copy.deepcopy(data)

    # b. age-correct
    moments.age_correct_data(par,data)
    
    # c. cov(dY_t, dY_t-1)         
    dnow = data.logY[2:,:] - data.logY[1:-1,:]
    dpre = data.logY[1:-1,:] - data.logY[0:-2,:]
    c = np.cov(dnow.reshape(-1) , dpre.reshape(-1))[0,1]
    model.sigma_xi = np.sqrt( -c )

    # d. cov(dY_t, y_t+1 - y_t-2) = cov(dY_t, dY_t+1 + dY_t + dY_t-1)
    dnow2 = data.logY[2:-1,:] - data.logY[1:-2,:] # y_t - y_t-1
    dbro = data.logY[3:,:] - data.logY[:-3,:]     # y_t+1 - y_t-2
    c = np.cov( dnow2.reshape(-1) , dbro.reshape(-1) )[0,1]
    model.sigma_psi = np.sqrt( c )

    print(f'sigma_psi = {model.sigma_psi:8.5f}, sigma_xi = {model.sigma_xi:8.5f}')
    
    # e. initial
    model.P_ini_std = np.sqrt(np.var(data.logY[0,:]) - model.sigma_xi**2 - model.sigma_psi**2)

    print(f'P_ini_std = {model.P_ini_std:8.4f}')

    return model

def simulate(par,model,data_in,seed=1917,rng_state=None):

    if not seed == None: 
        np.random.seed(seed)
    else:
        np.random.set_state(rng_state)
        
    # a. allocate
    data = LRT.dataStruct()
    data.logY = np.zeros((par.T,par.simN))
    data.logP = np.zeros((par.T,par.simN))

    # b. random
    psi = np.random.normal(-0.5*model.sigma_psi**2,model.sigma_psi,(par.T,par.simN))
    xi = np.random.normal(-0.5*model.sigma_xi**2,model.sigma_xi,(par.T,par.simN))
    logP0 = np.mean(data.logY[0,:]) + np.random.normal(-0.5*model.P_ini_std**2,model.P_ini_std,(1,par.simN)) 

    # c. simulate
    for t in range(par.T):    
        if t == 0:          
            data.logP[t,:] = logP0 + psi[t,:]
        else:
            Y = np.mean(np.exp(data_in.logY[t,:]))
            Y_lag = np.mean(np.exp(data_in.logY[t-1,:]))
            logGamma = np.log(Y/Y_lag)
            data.logP[t,:] = logGamma + data.logP[t-1,:] + psi[t,:]        
        data.logY[t,:] = data.logP[t,:] + xi[t,:]   

    return data