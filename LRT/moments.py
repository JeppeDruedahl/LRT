import copy
import time
import itertools as it 

import numpy as np
import pandas as pd
import scipy.stats

import matplotlib.pyplot as plt

from . import figfuns

class momentStruct:
    def __init__(self):
        pass

###############
# AGE CORRECT #
###############

def age_dummies(par,data):

    # a. initialize
    data.age_dummies = np.nan*np.ones(par.T)

    # b. calculate
    for t in range(par.T):

        # i. select
        I = np.isnan(data.logY[t,:]) == False
        
        # ii. calculate
        if np.sum(I) > 0:
            data.age_dummies[t] = np.mean(data.logY[t,I])
            
def age_correct_data(par,data):

    # a. age averages in logs
    age_dummies(par, data)

    # b. subtract averages
    for t in range(par.T):
        data.logY[t,:] -= data.age_dummies[t]

def age_decorrect_data(par,data):

    # a. subtract averages
    for t in range(par.T):
        data.logY[t,:] += data.age_dummies[t]

    # b. delete age dummies
    delattr(data,'age_dummies')


######################
# BASIC CALCULATIONS #
######################

def recent_earnings(par,data):

    # a. allocate
    data.RE = np.nan*np.ones(data.logY.shape)

    # b. time loop
    t0 = par.RE_min_obs # we need at least this many anyway
    for t in range(t0,par.T):

        # i. window
        tbeg = np.fmax(t-par.RE_lags,0)
        RE = np.zeros(data.logY.shape[1])
        num_obs = np.zeros(data.logY.shape[1])
        for tnow in range(tbeg,t):

            # include all not NaN
            I = np.isnan(data.logY[tnow,:]) == False
            RE[I] += data.Y[tnow,I] / np.exp(data.age_dummies[tnow])
            num_obs[I] += 1

        # ii. mean (if in base sample)
        I = (num_obs >= par.RE_min_obs)
        data.RE[t,I] = RE[I] /  num_obs[I]

def calc_dlogY(par,data,k):

    # a. initialize at NaN
    dlogY = np.nan*np.ones(data.logY.shape)

    # b. determine time period
    tmin = np.fmax(0,-k)
    tmax = np.fmin(par.T,par.T-k)

    # c. time loop
    for t in range(tmin,tmax):

        # i. now and future
        now = data.logY[t,:]
        future = data.logY[t+k,:]

        # ii. both incomes observed and non-zero
        I = (np.isnan(data.logY[t,:]) == False) & (np.isnan(data.logY[t+k,:]) == False)

        # ii. return
        dlogY[t,I] = future[I] - now[I]

    return dlogY

def calc_dlogRE(par,data,k):

    # a. initialize at NaN
    dlogRE = np.nan*np.ones(data.logY.shape)

    # b. determine time period
    tmin = np.fmax(0,-k)
    tmax = np.fmin(par.T,par.T-k)

    # c. time loop
    for t in range(tmin,tmax):

        # i. now and future
        now = data.RE[t,:]
        future = data.RE[t+k,:]

        # ii. both incomes observed and non-zero
        I = (np.isnan(now)==0) & (np.isnan(future)==0) 
        I[I] &= (now[I]>0.) & (future[I]>0.)

        # iii. return
        dlogRE[t,I] = np.log(future[I]) - np.log(now[I])

    return dlogRE


###############
# PERCENTILES #
###############

def perc_on_sorted(x,p):

    assert p >= 0
    assert p <= 1

    if len(x) == 0: 
        return np.nan

    # a. numerically safe
    if p <= 0:
        return x[0]
    elif p >= 1:
        return x[-1]

    # b. size of vector
    N = np.size(x)

    # c. neighbor indexes
    i1 = np.int(np.floor(p*(N-1)))

    if i1 == N-1:
        i1 = i1-1
        i2 = i1+1
    else:
        i2 = i1+1

    # d. neighbors percentiles
    p1 = i1/(N-1)
    p2 = i2/(N-1)

    # e. weight
    w = (p-p1)/(p2-p1)

    # f. retun
    return x[i1] + w*(x[i2]-x[i1])

def RE_percs(par,data,age_grp_str,percs_grp_str):

    # a. load
    age_grp = par.age_grps[age_grp_str]
    perc_grp = par.perc_grps[percs_grp_str]

    # b. initialize without a group
    RE_perc = -2*np.ones(data.logY.shape,dtype=np.int)

    # c. loop over age groups
    bounds = [0.,0.] # initialize
    
    # to avoid warnings
    RE = data.RE
    I = (np.isnan(RE) == 1)
    RE[I] = -np.inf
    
    for age_bounds in age_grp:

        # i. view of RE for current age group
        RE_now = RE[age_bounds[0]-par.agemin:age_bounds[1]-par.agemin+1,:]
        RE_perc_now = RE_perc[age_bounds[0]-par.agemin:age_bounds[1]-par.agemin+1,:]

        # ii. find non-nan
        I = (RE_now != -np.inf)

        # iii. set RE_perc to -1 for those with RE = 0
        J = (RE_now == 0)
        RE_perc_now[I&J] = -1

        # iv. non-zero RE percentile groups
        I[I] &= (RE_now[I] > 0)

        # v. sorted copy of RE to work with
        RE_sorted = np.sort(RE_now[I].ravel())

        for iperc, perc_bounds in enumerate(perc_grp):

            bounds[0] = perc_on_sorted(RE_sorted,perc_bounds[0]/100.)
            bounds[1] = perc_on_sorted(RE_sorted,perc_bounds[1]/100.)

            # note: this could cause a warning due to geg of nan
            J = (I == 1) & (RE_now >= bounds[0]) & (RE_now <= bounds[1])
            RE_perc_now[J] = iperc

    # d. output
    data.RE_perc[(age_grp_str, percs_grp_str)] = RE_perc


################
# MAIN MOMENTS #
################

def allocate_describe(dictnow,data,age_grp_str,perc_grp_str,varname_str,k,size):

    mean = np.nan*np.ones(size)
    var  = np.nan*np.ones(size)
    skew = np.nan*np.ones(size)
    kurt = np.nan*np.ones(size)

    dictnow[(varname_str, age_grp_str, perc_grp_str, 'mean', k)] = mean
    dictnow[(varname_str, age_grp_str, perc_grp_str, 'var' , k)] = var
    dictnow[(varname_str, age_grp_str, perc_grp_str, 'skew', k)] = skew
    dictnow[(varname_str, age_grp_str, perc_grp_str, 'kurt', k)] = kurt

    return mean, var, skew, kurt

def mainmoments(par,data,age_grp_str,perc_grp_str):

    # a. create dictionaries
    data.moments.levels = dict()
    data.moments.changes = dict()
    data.moments.autocorr = dict()
    data.moments.cov = dict()
    data.moments.cov_YdY = dict() 

    # b. load
    perc_grp = par.perc_grps[perc_grp_str]
    data_perc = data.RE_perc[(age_grp_str,perc_grp_str)]
    age_grp = par.age_grps[age_grp_str]

    # c. levels and autocorrelation, all ages
    dictnow = data.moments.levels
    varname_str = 'logY'
    size = (par.T)

    # covariances between logY at different ages
    C = np.cov(data.logY)
    for t in range(par.T): 
        data.moments.cov[('logY','cov',t)] = C[t:,t] # save lower diagonal only (time-space-continuum and all...)
        data.moments.cov[('logY','age',t)] = data.age[t:,0] # doesn't matter which column we use, so 0 

    mean, var, skew, kurt = allocate_describe(dictnow, data, 'all', 'off', varname_str, 0, size)
    for k in par.K_autocorr:
        autocorr = np.nan*np.ones(size)
        data.moments.autocorr[('dlogY', 'all','off', k)] = autocorr
        autocorr_level = np.nan*np.ones(size)
        data.moments.autocorr[('logY', 'all', 'off', k)] = autocorr_level

    for t in range(par.T):

        # i. pick data
        datavar = data.logY[t,:]
        I = (np.isnan(data.logY[t,:]) == False)
        datavar_masked = datavar[I]

        # ii. calculate moments
        mean[t] = np.mean(datavar_masked) + data.age_dummies[t]
        var[t]  = np.var(datavar_masked)
        skew[t] = scipy.stats.skew(datavar_masked)
        kurt[t] = scipy.stats.kurtosis(datavar_masked,fisher=False)

        # autocorrelation
        now = data.dlogY[1][t,:]
        now_level = data.logY[t,:]
        for k in par.K_autocorr:

            tfut = t+k
            if tfut >= par.T:
                continue

            # future
            fut = data.dlogY[1][tfut,:]
            fut_level = data.logY[tfut,:]

            # selected
            I = (np.isnan(now)==False) & (np.isnan(fut)==False)
            J = (np.isnan(now_level)==False) & (np.isnan(fut_level)==False)

            # calculate
            if np.any(I):
                autocorr = data.moments.autocorr[('dlogY', 'all','off', k)]
                autocorr[t] = np.corrcoef(now[I],fut[I])[1,0]
            if np.any(J): 
                autocorr = data.moments.autocorr[('logY', 'all','off', k)]
                autocorr[t] = np.corrcoef(now_level[J], fut_level[J])[1,0]

    # d. changes, all ages
    dictnow = data.moments.changes
    varname_str = 'dlogY'
    size = (par.T)
    k_list = par.K_leads
    
    for k in k_list:
            
        mean, var, skew, kurt = allocate_describe(dictnow, data, 'all', 'off', varname_str, k, size)
        for t in range(par.T):
            
            if t+k >= par.T:
                continue
            
            # i. pick data
            datavar = data.dlogY[k][t,:]
            I = np.isnan(datavar) == False
            datavar_masked = datavar[I]
    
            # ii. calculate moments
            mean[t] = np.mean(datavar_masked) + data.age_dummies[t+k] - data.age_dummies[t]
            var[t]  = np.var(datavar_masked)
            skew[t] = scipy.stats.skew(datavar_masked)
            kurt[t] = scipy.stats.kurtosis(datavar_masked,fisher=False)
    
    # e. leads and lags, age groups
    varname_str_list = par.mainmoments_varname_str_list 
    size = (len(perc_grp),len(age_grp))

    # allocate
    for varname_str in varname_str_list:

        if varname_str == 'logYk':
            k_list = par.K_leads_and_lags
            dictnow = data.moments.levels
        elif varname_str == 'dlogY':
            k_list = par.K_leads
            dictnow = data.moments.changes
        elif varname_str == 'dlogRE':
            k_list = par.K_leads_and_lags_long
            dictnow = data.moments.changes

        for k in k_list:
            allocate_describe(dictnow,data,age_grp_str,perc_grp_str,varname_str,k,size)

    varname_str = 'dlogY'
    for k in par.K_autocorr:
        autocorr = np.nan*np.ones(size)
        data.moments.autocorr[('dlogY', age_grp_str, perc_grp_str, k)] = autocorr

    # loop over age groups
    for iage, age_bounds in enumerate(age_grp):

        # i. slice
        data_perc_now = data_perc[age_bounds[0]-par.agemin:age_bounds[1]-par.agemin+1,:]
        data_logY_now = data.logY[age_bounds[0]-par.agemin:age_bounds[1]-par.agemin+1,:]
        J = np.isnan(data_logY_now) == False

        # ii. variable loop
        for varname_str in varname_str_list:

            if varname_str == 'logYk':
                k_list = par.K_leads_and_lags
                dictnow = data.moments.levels
            elif varname_str == 'dlogY':
                k_list = par.K_leads
                dictnow = data.moments.changes
            elif varname_str == 'dlogRE':
                k_list = par.K_leads_and_lags_long
                dictnow = data.moments.changes

            # lead and lag loop
            for k in k_list:

                # o. find output
                mean = dictnow[(varname_str, age_grp_str, perc_grp_str, 'mean', k)]
                var = dictnow[(varname_str, age_grp_str, perc_grp_str, 'var', k)]
                skew = dictnow[(varname_str, age_grp_str, perc_grp_str, 'skew', k)]
                kurt = dictnow[(varname_str, age_grp_str, perc_grp_str, 'kurt', k)]

                # oo. find data
                datavar_dict = getattr(data,varname_str)
                datavar = datavar_dict[k]
                datavar_now = datavar[age_bounds[0]-par.agemin:age_bounds[1]-par.agemin+1,:]

                # ooo. percentile groups
                for iperc in range(len(perc_grp)):

                    # select
                    I = data_perc_now == iperc
                    I &= J
                    I &= np.isnan(datavar_now) == 0
                    datavar_masked = np.ravel(datavar_now[I])

                    # calculate
                    if datavar_masked.size > 0:

                        mean[iperc,iage] = np.mean(datavar_masked) #+ data.age_dummies[t+k] - data.age_dummies[t]
                        var[iperc,iage]  = np.var(datavar_masked)
                        skew[iperc,iage] = scipy.stats.skew(datavar_masked)
                        kurt[iperc,iage] = scipy.stats.kurtosis(datavar_masked,
                                                                fisher=False)

        # iii. autocorrelation
        varname_str = 'dlogY'
        datavar_dict = getattr(data,varname_str)
        datavar = datavar_dict[1]

        for k in par.K_autocorr:

            # o. slice
            maxbound = np.fmin(age_bounds[1]-par.agemin,par.T-k-1)
            minbound = age_bounds[0]-par.agemin
            if minbound > maxbound:
                continue

            data_perc_now = data_perc[minbound:maxbound+1,:]
            
            # oo. now and future
            now = datavar[minbound:maxbound+1,:]
            fut = datavar[minbound+k:maxbound+k+1,:]
            
            # ooo. percentile groups
            for iperc in range(len(perc_grp)):

                I = data_perc_now == iperc
                I &= (np.isnan(now) == 0) & (np.isnan(fut) == 0)

                if np.any(I):
                    autocorr = data.moments.autocorr[('dlogY', age_grp_str, perc_grp_str, k)]
                    autocorr[iperc,iage] = np.corrcoef(now[I],fut[I])[1,0]

    # f. covariance of levels and future earnings growth 
    tt = par.cov_YdY_tt 
    ss = par.cov_YdY_ss
    kk = par.cov_YdY_kk
    dictnow = data.moments.cov_YdY

    for t,s in it.product(tt,ss): 
        if (t+s >= par.T): 
            continue

        yt = data.logY[t,:]
        yts = data.logY[t+s,:]

        dictnow[('cov_YdY',t,s)] = np.nan*np.zeros(shape=(len(kk)))
        for i,k in enumerate(kk): 
            if (t+s+k >= par.T): 
                continue 
            ytsk = data.logY[t+s+k,:]
            dictnow[('cov_YdY',t,s)][i] = np.cov(yt,ytsk-yts)[0,1]


################
# HETEROGENOUS #
################

def heterogenous(par, data):

    I = (data.RE[-1,:] > 0) & (data.RE[par.RE_min_obs,:] > 0)
    data.moments.heterogenous[('dlogY',0)] = np.log(data.RE[-1,I]) - np.log(data.RE[par.RE_min_obs,I])

    data.moments.heterogenous[('std_dlogY',1)] = np.nanstd(data.dlogY[1],axis=0)
    data.moments.heterogenous[('std_dlogRE',5)] = np.nanstd(data.dlogRE[5],axis=0)
    _T, N = data.logY.shape

    # autocorr of levels
    autocorr_level = np.ones(N)
    for i in range(N): 
        logY = data.logY[:,i]
        now = logY[1:-1]
        fut = logY[0:-2]
        autocorr_level[i] = np.corrcoef(now,fut)[1,0]
    data.moments.heterogenous[('autocorr_level',1)] = autocorr_level

    # autocorr of growth rates 
    autocorr = np.ones(N)
    for i in range(N):
        dlogY = data.dlogY[1][:,i]
        now = dlogY[1:-1]        
        fut = dlogY[0:-2]
        autocorr[i] = np.corrcoef(now,fut)[1,0]
    data.moments.heterogenous[('autocorr',1)] = autocorr

    
#######
# ALL #
#######

def calc_all(par,data,printprogress=False):

    if printprogress:
        print('calculating all moments:')
    data.moments = momentStruct()

    ###############
    # age correct #
    ###############
    
    # a. construct age
    if hasattr(data,'age') == False:
        data.age = np.ones(data.logY.shape)
        for t in range(par.T):
            data.age[t,:] = par.agemin + t

    # b. calculate Y
    if hasattr(data,'Y') == False:        
        data.Y = np.nan*np.ones(data.logY.shape) # missings: for m==1
        I = np.isnan(data.logY) == False      
        data.Y[I] = np.exp(data.logY[I]) 
    
    # c. correct earnings for age effects
    if hasattr(data,'age_dummies') == False:
        age_correct_data(par,data)
    

    ######################
    # basic computations #
    ######################

    t1 = time.time()
    recent_earnings(par, data)

    if printprogress:
        print(f' - basic computations ({time.time()-t1:3.1f} secs)')

    t1 = time.time()
    data.dlogY = dict()
    data.dlogRE = dict()
    for k in par.K_leads:
        data.dlogY[k] = calc_dlogY(par,data,k)
    for k in par.K_dlogRE_leads:
        data.dlogRE[k] = calc_dlogRE(par,data,k)
    if printprogress:
        print(f' - recent earnings ({time.time()-t1:3.1f} secs)')


    #############################################
    # recent earnings percentiles by age groups #
    #############################################

    t1 = time.time()
    data.RE_perc = dict()
    RE_percs(par, data, 'detail', 'RE')
    if printprogress:    
        print(f' - recent earnings percentiles ({time.time()-t1:3.1f} secs)')


    ###############################
    # moments of logY, dlogY, dlogRE #
    ###############################

    t1 = time.time()
    mainmoments(par, data, age_grp_str='detail', perc_grp_str='RE')
    if printprogress:     
        print(f' - main-moments ({time.time()-t1:3.1f} secs)')

    ####################################
    # individual heterogeneity in risk #
    ####################################

    t1 = time.time()    
    data.moments.heterogenous = dict()
    heterogenous(par, data)
    if printprogress:         
        print(f' - heterogenous ({time.time()-t1:3.1f} secs)')

    # over and out
    age_decorrect_data(par,data)


####################
# ROBUST MOMEMENTS #
#################### 

def calc_robust_moment(y, mom, diff): 
    
    if diff: 
        assert y.shape[0] == 2
        vec = y[1,:] - y[0,:]
    else: 
        vec = y[0,:]
    I = np.isnan(vec) == False

    vec = vec[I]
    if mom=='mean': 
        ans = np.mean(vec)
    elif mom=='var': 
        pp = np.percentile(vec,q=[10.,90.])
        ans = pp[1] - pp[0]
    elif mom=='skew': 
        pp = np.percentile(vec[I], q=[10.,50.,90.])
        ans = ((pp[2]-pp[1])-(pp[1]-pp[0]))/(pp[2]-pp[0])
    elif mom=='kurt': 
        pp = np.percentile(vec, q=[12.5,25.,37.5,62.5,75.,87.5])
        ans = ((pp[5]-pp[3]) + (pp[2]-pp[0])) / (pp[4]-pp[1])
            
    return ans

def calc_robust_moments_all_ages(dat,mom,diff,Tup): 
    return [calc_robust_moment(dat[t:t+2,:],mom,diff) for t in range(Tup)]


##########
# WITHIN #
##########

def weight_group_vec_by_counts(model,group_vec_in,t): 
    
    # unpack
    N = model.G[t].shape[0]
    counts = pd.DataFrame(model.G[t])[0].value_counts().sort_index().values
    group_vec = copy.deepcopy(group_vec_in)

    # handle missings 
    for i in range(len(group_vec)): 
        if np.isnan(group_vec[i]): 
            counts[i] = 0.
            group_vec[i] = 0.

    # sum
    N = np.sum(counts)            

    return np.dot(counts, group_vec)/N

def analytic_jth_moment_of_kyr_growth(model,t,k,j): 
    
    # a. income now/then 
    ynow  = model.ypred_G[t]
    yplus = model.ypred_G[t+k]

    # b. transition matrix 
    T = model.trans[t+1]
    for this_k in np.arange(2, k+1): 
        T = np.dot(T, model.trans[t+this_k])

    diff = -np.subtract.outer(ynow,yplus)
    mean_diff = np.sum(diff * T, axis=1)
    diff_demeaned = ((diff.transpose()-mean_diff).transpose()) 

    # c. variance
    y = np.sum((diff_demeaned**2) * T, axis=1)

    # d. skew/kurt
    if j >= 3: 
        s = np.sqrt(y)
        Z = (diff_demeaned.transpose()/s).transpose()
        integrand = Z**j
        y = np.sum(integrand * T, axis=1)
    
    return y 

def compute_jth_central_moment_of_vec(vec,moment):

    if moment == 2: 
        y = vec.var()
    elif moment == 3: 
        y = scipy.stats.skew(vec)
    elif moment == 4: 
        y = scipy.stats.kurtosis(vec)

    return y