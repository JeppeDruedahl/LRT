import numpy as np 
import copy 
import os 

from . import LRT

def count_num_obs_cens(T_row, I_row): 
    if I_row.any(): 
        return T_row[I_row].sum()
    else: 
        return 0
    
import os 

def write_model_to_disk(out_dir, model, min_obs_per_cell=5, DOCENSOR=True): 

    # common fdormat for printing floats 
    float_fmt = '%16.11f'
    
    if not os.path.isdir(out_dir): 
        print(f'Directory, {out_dir}, not found. Creating it')
        os.mkdir(out_dir)

    # a. number of leafs
    np.savetxt(out_dir + "/num_leafs.txt",model.num_leafs,delimiter=',',fmt='%d')

    # b. initial grouping
    prob_G_ini = np.ones(model.num_leafs[0])
    N = model.G[0].shape[0] # number of individuals in dataset 
    for i in range(model.num_leafs[0]):
        # i. fill out 
        I = (model.G[0] == model.uniqueG[0][i])
        prob_G_ini[i] = I.sum()/N

        # ii. assertion 
        assert np.nanmin(prob_G_ini*N) > min_obs_per_cell, f'Marginals are not anonymous!!'

    np.savetxt(out_dir + "/prob_G_ini.txt",prob_G_ini,delimiter=',',fmt=float_fmt)

    # c. save predicted income in a single file 
    T = len(model.ypred_G) # number of ages 
    assert len(model.G) == T , 'Unexpected disagreement of dimensions between G and ypred_G'
    maxNG = np.max(model.num_leafs)
    ypred_G = np.nan * np.empty((maxNG, T))
    for t in range(T): 
        ypred_G[:model.num_leafs[t], t] = model.ypred_G[t]
    filename = f'{out_dir}/ypred_G.txt'
    np.savetxt(filename,ypred_G,delimiter=',',fmt=float_fmt)    

    # d. censor transition matrix 
    if DOCENSOR: 
        T_c, T_obs_c = censor_trans_mat(model)
        assert all_transitions_with_anonymous_obs(T_obs_c, min_obs_per_cell), f'T_obs_c contains non-anonymous cells'
    else: 
        T_c = model.trans
        T_obs_c = model.trans_obs

    # e. save censored transition matrices
    for i,this_T in enumerate(T_c):  
        if this_T == []: 
            continue
            
        if DOCENSOR: 
            t = i+1
        else: 
            t = i 
        
        trans = np.copy(this_T)
        filename = f"{out_dir}/trans_t{t:d}.txt"
        np.savetxt(filename,trans,delimiter=',',fmt=float_fmt) 

def find_second_lowest_col_not_already_censored(T_row, I_already_being_censored, PRINTDETAILS=False): 
    """
        Returns the index of (one of) the column that is not already being censored.

        ARGS:
            T_row: row-vector of observation counts for cells in the transition matrix
            I_already_being_censored: list of logicals indicating rows that are already being censored. 
        RETURNS:
            (scalar, int) Index for (one of) the column(s) containing the second-lowest number of obs. 
            and which was not already being censored (as indicated by I_already_being_censored). 
    """
    ii = np.arange(T_row.shape[0])
    
    # remove the cells that we should disregard for either 
    #   - containing zero obs. 
    #   - already being censored 
    I_also_not_these = T_row == 0
    I = (I_already_being_censored == False) & (I_also_not_these == False)
    T_short = T_row[I]
    ii = ii[I]
    
    # find smallest obs. count amont the remainder 
    I_smallest = (T_short == np.min(T_short))
    ii_smallest = ii[I_smallest] # column number of that cell 
    
    # if there are multiple, pick one at random 
    if len(ii_smallest) == 1: 
        i = ii_smallest[0]
    else: 
        i = np.random.choice(ii_smallest)

    
    if PRINTDETAILS:
        print(f'There are {ii_smallest.shape[0]} to choose from; returning i={i} (of {ii_smallest})')
        
    return i

def all_transitions_with_anonymous_obs(T_obs, min_obs_per_cell=5): 
    """
        Verify that all transition matrix cells are properly anonymized

        INPUT: matrix with number of obs. in each transition 
        RETURNS: bool, True if no cells have 0 < obs <= min_obs_per_cell
    """
    ret = True
    
    for t in range(len(T_obs)): 

        T_ = T_obs[t]

        # 1. create I: indicator for (T_>0) & (np.isnan(T_) == False)

        # 1.a initialize to False (willk cover all np.isnan)
        I = np.zeros(T_.shape, dtype='bool') # default=False: do not select nan cells

        # 1.b set to True where T_>0 (for cells that are not nan)
        I[np.isnan(T_) == False] = (T_[np.isnan(T_) == False] > 0)

        # 2. make assertions 
        if not np.any(I):
            print(f'Weird, at t={t}, no cell has obs > 0?')
            ret = False 
            
        lowest_obs_cell = np.min(T_[I])
        if not (lowest_obs_cell > min_obs_per_cell): 
            print(f'The censored cell matrix has a cell with {lowest_obs_cell} observations! Not anonymized')
            ret = False
            
    return ret

def censor_trans_mat(model, PRINTDETAILS=False): 
    """
        Censors the transitions so that no cell has 0<obs<=5. Furthermore, if any row would have 5 or fewer 
        obs. censored (so that exact numbers can be computed from marginal counts), we censor the cell with
        the second-lowest observation count (or one of those chosen at random in case there are multiple 
        cells with exactly the same number of observations). 
    
        INPUT: 
            model: must have "trans" and "trans_mat"
            
        RETURNS: 2 lists of: 
            Transition matrix that has been censored
            Obs. count matrix with same censoring scheme
    
    """
    
    if PRINTDETAILS: 
        print(f'Censoring model: {model.name}.')

    T_censored = []
    T_obs_censored = [] 
    
    total_cells_censored = [] 
    total_cells = [] 
    total_obs_censored = []
    total_obs = []

    for t in np.arange(1, len(model.trans)): 
        Tm = model.trans_obs[t]

        this_T_censored     = copy.deepcopy(model.trans[t])
        this_T_obs_censored = copy.deepcopy(Tm)

        # initialize indicator for censoring 
        I = (Tm>0) & (Tm<=5)
        num_cells_cens_this_row = np.sum(I, 1)

        for i_row in range(I.shape[0]): 
            # point to current row 
            this_row_obs = Tm[i_row,:]

            # how many obs. are about to be censored? 
            num_obs_cens_this_row = count_num_obs_cens(Tm[i_row,:], I[i_row,:])        

            # if the default censoring scheme has not lead to enough obs. being censored to ensure anonymity 
            if (num_cells_cens_this_row[i_row] == 1) or ((num_obs_cens_this_row <= 5) and (num_obs_cens_this_row > 0)): 
                
                # also censor (one of the) second-lowest rows
                # (i.e. i_also_cens points to a column that would not previously have been censored)
                i_also_cens = find_second_lowest_col_not_already_censored(Tm[i_row,:], I[i_row,:], PRINTDETAILS)

                assert np.sum(Tm[i_row, i_also_cens]) > 5, f'Unexpectedly low number of cells added for additional censoring, this will not yield anonymity.'
                assert not I[i_row, i_also_cens] , f'Found a cell to censor which has allready been assigned for censoring. Unexpected!'

                if PRINTDETAILS:
                    print(f'Updating cell ({i_row}, {i_also_cens})')
                I[i_row, i_also_cens] = True 

            else: 
                if PRINTDETAILS: 
                    print(f'Nothing further to censor in row {i_row}')


        # verify that we never censor {1,2,3,4,5} obs. 
        num_being_censored = np.sum(Tm*I, 1)
        assert not np.any((num_being_censored>0) & (num_being_censored<=5)), f'Too few obs. being censored.'
        
        # update 
        total_cells_censored.append(I.sum())
        total_cells.append(np.prod(I.shape))
        total_obs_censored.append(num_being_censored.sum())
        total_obs.append(np.sum(Tm))

        # do the censoring 
        this_T_censored[I] = np.nan
        this_T_obs_censored[I] = np.nan 

        if PRINTDETAILS: 
            print(f't={t:2d}: #censored cells: {np.sum(I, 1)}') 

        T_censored.append(this_T_censored)
        T_obs_censored.append(this_T_obs_censored)
        
    
    tc = np.sum(total_obs_censored)
    to = np.sum(total_obs) 
    print(f'Total censoring: {tc/to* 100.0: 5.4f}% ({tc:7.0f}/{to:7.0f})')
    tc = np.sum(total_cells_censored)
    to = np.sum(total_cells)
    print(f'Cells censored: {tc/to * 100.0: 5.4f}% ({tc:7.0f}/{to:7.0f})')
    
    return T_censored, T_obs_censored

def read_censored_data(in_dir, T=30, N=292840, marker='o', color='#1f77b4'): 
    """
        Read the matrices expected to be found in the "in_dir" directory: 
            num_leafs.txt
            prob_G_ini.txt: distribution over discrete states at t=0
            ypred_G_t{X}.txt: mean incomes for each age 
            trans_t{X}.txt: transition matrix between age X-1 and X
    """

    assert os.path.isdir(in_dir), f'Path, "{in_dir}", is not a directory.'

    m_load = LRT.modelStruct()
    m_load.name = 'LRT'
    m_load.marker = marker
    m_load.color = color

    # number of discrete groups at each age category 
    m_load.num_leafs = np.loadtxt(f'{in_dir}/num_leafs.txt', delimiter=',')
    assert np.all(np.isclose(m_load.num_leafs, np.round(m_load.num_leafs))), f'num_leafs are not all integers'
    m_load.num_leafs = m_load.num_leafs.astype(int)

    # initial distribution
    m_load.prob_G_ini = np.loadtxt(f'{in_dir}/prob_G_ini.txt', delimiter=',')
    fname = f'{in_dir}/ypred_G.txt'
    ypred_G = np.loadtxt(fname, delimiter=',')

    m_load.ypred_G = [] 
    for t in range(T): 
        m_load.ypred_G.append(ypred_G[:m_load.num_leafs[t], t])
        assert not np.any(np.isnan(m_load.ypred_G[t])), f'{t}: Found NaNs in predicted incomes'
        if m_load.num_leafs[t] < ypred_G.shape[0]: # last rows should be NaNs
            remainder = ypred_G[m_load.num_leafs[t]:, t]
            assert np.all(np.isnan(remainder)), f'{t}: Found not NaNs where NaNs were expected (rows after num_leafs={m_load.num_leafs[t]})'

    # transition matrices 
    m_load.trans = [] 
    m_load.trans.append([]) # t=0 is empty
    for t in range(1, T): 
        # load transition matrix from t-1 to t
        fname = f'{in_dir}/trans_t{t:d}.txt'
        Trans = np.loadtxt(fname, delimiter=',')

        # handle nans
        Trans_normalized = renormalize_censored_transition_matrix(Trans)
        
        # put in list 
        m_load.trans.append(Trans_normalized)

        # verify sensibility
        assert np.all(m_load.trans[t].shape == (m_load.num_leafs[t-1], m_load.num_leafs[t])), f'Trans({t-1}, {t}): unexpected dimensions'



    return m_load


def renormalize_censored_transition_matrix(T_in): 
    """
        For all rows in the transition matrix, replace the nans with an equal share of the 
        probability mass that is missing for that row. 
    """
    Tmat = copy.deepcopy(T_in)
    for i_row in range(Tmat.shape[0]): 
        I = np.isnan(Tmat[i_row, :])
        if np.any(I): 
            n_cells = np.sum(I)
            missing_mass = 1.0 - np.sum(Tmat[i_row, I == False])
            Tmat[i_row, I] = missing_mass / n_cells
    
    return Tmat 
    

def simulate(model,N=292840,T=30,seed=1917,rng_state=None): 
    
    if not seed == None: 
        np.random.seed(seed)
    else:
        np.random.set_state(rng_state)

    # a. seutp
    data = LRT.dataStruct() 
    data.G = np.empty((T,N), dtype=np.int)
    data.logY = np.empty((T,N))
    
    # b. simulate
    for t in range(T): 
        if t == 0: 
            
            # must be in [0,1,...,model.num_leafs[0]-1] and not in model.uniqueG[0]
            data.G[0,:] = draw_initial_distribution_idx(model, N) #np.random.choice(model.num_leafs[0], N, p=model.prob_G_ini)
        
        else: 

            for i_past in range(model.num_leafs[t-1]): 

                I = data.G[t-1,:] == i_past
                num = np.sum(I)
                
                if num > 0:

                    pvec = model.trans[t][i_past,:]
                    data.G[t,I] = np.random.choice(model.num_leafs[t], num, p=pvec)
        
        data.logY[t,:] = model.ypred_G[t][data.G[t,:]]
    
    return data 

def draw_initial_distribution_idx(model, N, seed=1917, DETERMINISTICDRAWS=True): 
    """
        Draws an initial distribution (N obs.) of {0,1,...,num_leafs[0]-1} 
        using the probabilities in model.prob_G_ini. 
        
        ARGS: 
            model: must have members num_leafs and prob_G_ini
            N: number of obs. to draw
            seed: for the RNG 
            DETERMINISTICDRAWS: whether to just use prob_G_ini*N
            
        RETURNS: 
            G0_dist: N-vector of integers {0, 1, ..., model.num_leafs[0]-1}
    """
    
    if DETERMINISTICDRAWS: 
        # raw simulation, in floating points 
        raw_sim = model.prob_G_ini*N
        
        # exact counts as int
        counts = np.round(raw_sim).astype(int)
        assert np.sum(counts)==N, f'Rounding off the raw probabilities implied by prob_G_ini*N does not give N={N} observations. '

        # repeat each discrete type
        G0_dist = np.repeat(np.arange(model.num_leafs[0]), counts)

        # randomize the order 
        np.random.seed(seed)
        np.random.shuffle(G0_dist)
    else: 
        
        np.random.seed(seed)
        G0_dist = np.random.choice(model.num_leafs[0], N, model.prob_G_ini)
    
    return G0_dist 
