import os
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

color_list  = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf']
marker_list = ['o','x','d','v','^','<','>']
linestyle_list = ['--','-.',':',(0,(1,10)),(0,(5,10))]
plt.rc('axes', prop_cycle=(cycler('color',color_list)))
plt.rcParams.update({'font.size': 12,'font.family': 'STIXGeneral','mathtext.fontset': 'stix'})

gray = (65.0/255.0,68.0/255.0,81.0/255.0)

from . import moments

def new():

    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)
    ax = fig.add_subplot(1,1,1)

    ax.spines["top"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.spines["right"].set_linewidth(0.5)
    ax.spines["left"].set_linewidth(0.5)
    ax.tick_params(which='both', color=gray)
    ax.grid(True, zorder=0, color=gray)

    return fig, ax

def save(fig, ax, name):

    [line.set_zorder(3) for line in ax.lines]
    title = ax.get_title()
    ax.set_title('')
    fig.savefig(os.getcwd() + '/output/' + name + '.pdf')
    ax.set_title(title)

def ylabel_strings(mom_str,varname_str,k): 

    if varname_str == 'logY':
        ylabel_var_str = '$y_{it}$' 
    else:
        ylabel_var_str = '$y_{it+' + str(k) + '} - y_{it}$' 
    
    if mom_str=='mean': 
        moment_name_str = 'mean'
    elif mom_str=='var': 
        moment_name_str = 'variance'
    elif mom_str=='skew': 
        moment_name_str = 'skewness'
    elif mom_str=='kurt': 
        moment_name_str = 'kurtosis'
    elif mom_str=='autocorr': 
        moment_name_str = 'autocorrelation'

    return ylabel_var_str,moment_name_str

def age_profile(par,data,models,k,varname_str,prefix='',save_to_disk=True):

    if varname_str == 'logY':
        dictnow = data.moments.levels 
    else:
        dictnow = data.moments.changes

    for mom_str in ['mean','var','skew','kurt']:

        fig, ax = new()

        # plot
        x = np.arange(par.agemin,par.agemax+1)
        y = dictnow[(varname_str, 'all', 'off', mom_str, k)]
        ax.plot(x,y,'-s',color='black',label='data')
        
        for model in models: 
            if varname_str == 'logY':
                dictmodel = model.data.moments.levels 
            else:
                dictmodel = model.data.moments.changes

            y = dictmodel[(varname_str, 'all', 'off', mom_str, k)]
            ax.plot(x,y,'-',marker=model.marker,color=model.color,label=model.name)

        # details
        ax.set_xlabel('age')

        ylabel_var_str,moment_name_str = ylabel_strings(mom_str,varname_str,k)
        ax.set_ylabel(moment_name_str + ' of ' + ylabel_var_str)
        
        # title
        tit = ''
        if prefix != '':
            tit += prefix + '_'
        tit += varname_str + '_' + mom_str 
        if varname_str != 'logY':        
            tit += '_k' + str(k)
        ax.set_title(tit)

        # only legend on the mean graph 
        if mom_str == 'mean': 
            legend = ax.legend(loc='best',fontsize=10)
            legend.get_frame().set_edgecolor(gray)
            legend.get_frame().set_linewidth(.5)
                        
        # save
        if prefix != '' and save_to_disk:
            fig.tight_layout() 
            save(fig, ax, tit)
            plt.close()
        else: 
            plt.show()

def cov(par,data,model,modelname='model',agestep=10,prefix='',save_to_disk=True):
    
    # plot
    figdata,axdata = new() 
    figmodel,axmodel = new() 

    for t in range(par.T): 
        x = data.moments.cov[('logY','age',t)]
        y = data.moments.cov[('logY','cov',t)]
        axdata.plot(x,y)

    for t in range(par.T): 
        x = model.data.moments.cov[('logY','age',t)]
        y = model.data.moments.cov[('logY','cov',t)]
        axmodel.plot(x,y)

    axdata.set_xlabel('age')
    axdata.set_ylabel('covariance')   
    axmodel.set_xlabel('age')
    axmodel.set_ylabel('covariance')

    # common y axis limits 
    ydL,ydH = axdata.get_ylim()
    ymL,ymH = axmodel.get_ylim()
    ylim = [min([ydL,ymL]), max([ydH,ymH])]
    axdata.set_ylim(ylim)
    axmodel.set_ylim(ylim)

    # title
    tit = ''
    if prefix != '':
        tit += prefix + '_'
    tit_data = tit + 'cov_by_age_data'
    tit_model = tit + 'cov_by_age_model'

    axdata.set_title(tit_data)
    axmodel.set_title(tit_model)

    if prefix != '' and save_to_disk: 
        save(figdata, axdata, tit_data)
        save(figmodel, axmodel, tit_model)
        plt.close(figmodel)
        plt.close(figdata)
    else: 
        plt.show()
    
def autocorr(par,data,model,k_list,varname_str,modelname='model',prefix='',save_to_disk=True):

    # plot
    figdata, axdata = new()
    figmodel, axmodel = new()
    
    x = np.arange(par.agemin,par.agemax+1)
    for i,k in enumerate(k_list):
        ydata = data.moments.autocorr[(varname_str, 'all', 'off', k)]
        axdata.plot(x,ydata,marker=marker_list[i],color=color_list[i],label=f'order {k}')
        ymodel = model.data.moments.autocorr[(varname_str, 'all', 'off', k)]
        axmodel.plot(x,ymodel,marker=marker_list[i],color=color_list[i],label=f'order {k}')
    
    # title
    tit = ''
    if prefix !='':
        tit += prefix + '_'
    tit += 'autocorr_' + varname_str

    # limits
    y1_min, y1_max = axdata.get_ylim()
    y2_min, y2_max = axmodel.get_ylim()
    y_min = min([y1_min, y2_min])
    y_max = max([y1_max, y2_max])
    
    # detail
    for fig, ax, name in zip([figdata,figmodel],[axdata, axmodel],['data',modelname]):
            
        ax.set_title(name + '_' + tit)
        ax.set_ylim([y_min, y_max])    
        ax.set_xlabel('age')

        if name == 'data':
            legend = ax.legend(loc='best',fontsize=10)
            legend.get_frame().set_edgecolor(gray)
            legend.get_frame().set_linewidth(.5)
                    
        ylabel_var_str,moment_name_str = ylabel_strings('autocorr',varname_str,1)
        ax.set_ylabel(moment_name_str + ' of ' + ylabel_var_str)
        
        # save
        if save_to_disk: 
            fig.tight_layout() 
            save(fig,ax,name + '_' + tit)
            plt.close()

    if save_to_disk == False: 
        plt.show()
    
def cov_YdY(par,data,model,prefix='',save_to_disk=True): 
    
    figdata,axdata = new() 
    figmodel,axmodel = new() 

    for moments,ax,fig in zip([data.moments, model.data.moments],[axdata,axmodel],[figdata,figmodel]): 
        x = par.cov_YdY_kk

        i = 0
        for t in par.cov_YdY_tt: 
            for s in par.cov_YdY_ss: 
                if s+t >= par.T: 
                    continue 
                y = moments.cov_YdY[('cov_YdY',t,s)]
                ax.plot(x,y,marker=marker_list[i],color=color_list[i],label=f'$t = {t}, s = {s}$')
                i += 1

        ax.set_xlabel('$k$')
        ax.set_ylabel('cov$(y_t, y_{t+s+k} - y_{t+s})$')
        
        # avoid y labels hitting the left edge 
        fig.tight_layout() 

    legend = axdata.legend(loc='best',fontsize=10)
    legend.get_frame().set_edgecolor(gray)
    legend.get_frame().set_linewidth(.5)

    # set common y limits
    y1_min, y1_max = axdata.get_ylim()
    y2_min, y2_max = axmodel.get_ylim()
    y_min = min([y1_min, y2_min])
    y_max = max([y1_max, y2_max])

    for ax in [axdata,axmodel]: 
        ax.set_ylim([y_min,y_max])

    # title
    tit = ''
    if prefix != '': 
        tit += prefix + '_cov_YdY'
    tit_data = tit + '_data'
    tit_model = tit + '_model' 

    axdata.set_title(tit_data)
    axmodel.set_title(tit_model)

    if prefix != '' and save_to_disk: 
        save(figdata, axdata, tit_data)
        save(figmodel, axmodel, tit_model)
        plt.close(figdata)
        plt.close(figmodel)
    else: 
        plt.show()

def heterogenous(par,data,models,mom_str,k,xlabel,bounds,prefix='',save_to_disk=True):

    # plot
    fig, ax = new()   

    x = data.moments.heterogenous[(mom_str,k)]
    x = x[~np.isnan(x)]
    y = np.linspace(0,1,x.size)
    x = np.sort(x)
    ax.plot(x,y,'-',color='black',label='data')

    for i,model in enumerate(models): 

        x = model.data.moments.heterogenous[(mom_str,k)]
        x = x[~np.isnan(x)]
        y = np.linspace(0,1,x.size)
        x = np.sort(x)

        ax.plot(x,y,'-',color=model.color,linestyle=linestyle_list[i],label=model.name)

    # title
    tit = ''
    if prefix != '':
        tit += prefix + '_'
    tit += mom_str + f'_k{k}' 
    ax.set_title(tit)

    ax.set_xlabel(xlabel)
    ax.set_ylabel('')
    ax.set_xlim(bounds)

    # only legend on the dlogY graph 
    if mom_str == 'dlogY': 
        legend = ax.legend(loc='best',fontsize=10)
        legend.get_frame().set_edgecolor(gray)
        legend.get_frame().set_linewidth(.5)
                
    # save
    if prefix != '' and save_to_disk:
        fig.tight_layout()         
        save(fig, ax, tit)
        plt.close(fig)
    else: 
        plt.show()        

def REpercs_profile(par,data,model,k,mom_str,varname_str,age_grp_str,perc_grp_str,modelname='model',prefix='',save_to_disk=True):

    # plot
    figdata, axdata = new()
    dictdata = data.moments.changes
    
    figmodel, axmodel = new()
    dictmodel = model.data.moments.changes
    
    x = par.perc_x[perc_grp_str]
    ydata = dictdata[(varname_str, age_grp_str, perc_grp_str, mom_str, k)]
    ymodel = dictmodel[(varname_str, age_grp_str, perc_grp_str, mom_str, k)] 
    
    i = 0
    for iage, age_bounds in enumerate(par.age_grps[age_grp_str]):
        if np.mean(np.isnan(ydata[:,iage]))<1.0: 
            label = f'age {age_bounds[0]}-{age_bounds[1]}'
            axdata.plot(x,ydata[:,iage],marker=marker_list[i],color=color_list[i],label=label)
            axmodel.plot(x,ymodel[:,iage],marker=marker_list[i],color=color_list[i],label=label)
            i += 1

    # title
    tit = ''
    if prefix != '':
        tit += prefix + '_'
    if k < 0:
        tit += varname_str + '_' + mom_str + '_k_neg_' + str(-k)
    else:
        tit += varname_str + '_' + mom_str + '_k' + str(k)
    tit += '_age_' + age_grp_str
    
    # limits
    y1_min, y1_max = axdata.get_ylim()
    y2_min, y2_max = axmodel.get_ylim()
    y_min = min([y1_min, y2_min])
    y_max = max([y1_max, y2_max])

    # detail
    for fig, ax, name in zip([figdata,figmodel],[axdata, axmodel],['data',modelname]):
        
        ax.set_title(name + '_' + tit)
        ax.set_ylim([y_min, y_max])    
        ax.set_xlabel('percentiles of recent income')

        ylabel_var_str,moment_name_str = ylabel_strings(mom_str,varname_str,k)
        ax.set_ylabel(moment_name_str + ' of ' + ylabel_var_str)

        # only show legend on graph for mean 
        if mom_str == 'mean' and name == 'data': 
            legend = ax.legend(loc='best',fontsize=10)
            legend.get_frame().set_edgecolor(gray)
            legend.get_frame().set_linewidth(.5)
                        
        # save 
        if save_to_disk:
            fig.tight_layout() 
            save(fig, ax, name + '_' + tit)
            plt.close()
        else:
            plt.show()

def REpercs_profile_compare_models(par,data,models,k,mom_str,varname_str,age_grp_str,perc_grp_str,iage,prefix='',save_to_disk=True):

    age_bounds = par.age_grps[age_grp_str][iage]

    # plot
    fig, ax = new()
        
    x = par.perc_x[perc_grp_str]
    dictdata = data.moments.changes
    ydata = dictdata[(varname_str,age_grp_str,perc_grp_str,mom_str,k)]
    
    ax.plot(x, ydata[:,iage],'-s',label='data',color='black')

    for model in models: 

        dictmodel = model.data.moments.changes
        ymodel = dictmodel[(varname_str,age_grp_str,perc_grp_str,mom_str,k)] 
        ax.plot(x,ymodel[:,iage],'-',marker=model.marker,label=model.name,color=model.color)
        
    # title
    tit = ''
    if prefix != '':
        tit += prefix + '_'
    tit += f'age_{age_bounds[0]}_{age_bounds[1]}_'
    if k < 0:
        tit += varname_str + '_' + mom_str + '_k_neg_' + str(-k)
    else:
        tit += varname_str + '_' + mom_str + '_k' + str(k)
    tit += '_age_' + age_grp_str
    ax.set_title(tit)
    
    # details  
    ax.set_xlabel('percentiles of recent income')

    if iage == 0:
        legend = ax.legend(loc='best',fontsize=10)
        legend.get_frame().set_edgecolor(gray)
        legend.get_frame().set_linewidth(.5)
                
    ylabel_var_str,moment_name_str = ylabel_strings(mom_str,varname_str,k)
    ax.set_ylabel(moment_name_str + ' of ' + ylabel_var_str)

    # save
    if save_to_disk: 
        fig.tight_layout()
        save(fig,ax,tit)
        plt.close(fig)
    else: 
        plt.show()

def robust_age_profile(data_levels,models,mom,diff,prefix='',save_to_disk=False): 
    
    par = models[0].par 

    momNames = dict()
    momNames['mean'] = 'mean'
    momNames['var'] = 'P90-P10'
    momNames['skew'] = 'Kelley skewness'
    momNames['kurt'] = "Moors kurtosis"

    # changes or levels
    if diff: 
        yvar_lab = '$Y_{it+1} - Y_{it}$'
        yvar = 'dY'
        Tuppper = par.T-1
        ages = np.arange(par.agemin,par.agemax)
    else: 
        yvar_lab = '$Y_{it}$'
        yvar = 'Y'
        Tuppper = par.T
        ages = np.arange(par.agemin,par.agemax+1)

    # plot
    fig, ax = new()
    
    y = moments.calc_robust_moments_all_ages(data_levels.Y,mom,diff,Tuppper)
    plt.plot(ages,y,'-ks',label='data')
    
    for model in models: 
        y = moments.calc_robust_moments_all_ages(model.data.Y,mom,diff,Tuppper)
        plt.plot(ages,y,color=model.color,marker=model.marker,label=model.name)

    # details
    ax.legend()
    ax.set_xlabel('age')
    ax.set_ylabel(f'{momNames[mom]} of {yvar_lab}')

    # title
    tit = ''
    if prefix != '': 
        tit += prefix + '_' + yvar + '_' + mom

    if prefix != '' and save_to_disk: 
        fig.tight_layout()
        save(fig,ax,tit)
        plt.close(fig)
    else: 
        plt.show()

    ax.set_title(tit)

def moments_within_groups(data,model,marker_list,color_list,central_moment=2,kk=[1,5,10],prefix='',save_to_disc=True): 
    
    # unpack
    par = model.par
    Nk = len(kk)

    if central_moment == 2: 
        mom_name = 'var'
    elif central_moment == 3: 
        mom_name = 'skew'
    elif central_moment == 4: 
        mom_name = 'kurtosis'

    # compute
    y_uncond = np.empty((par.T-1,Nk))
    y_cond = np.empty((par.T-1,Nk)) 

    for ik,k in enumerate(kk):
        for t in range(par.T-k): 
            _y_cond = moments.analytic_jth_moment_of_kyr_growth(model,t,k,central_moment)
            y_cond[t,ik] = moments.weight_group_vec_by_counts(model,_y_cond,t)
            y_uncond[t,ik] = moments.compute_jth_central_moment_of_vec(model.data.logY[t+k,:]-model.data.logY[t,:],moment=central_moment)

    # plot
    for ik,k in enumerate(kk): 
        
        Tmax = par.T-k
        ages = np.arange(par.agemin, par.agemax)
        ages = ages[:Tmax]
        
        # plot
        fig,ax = new()
        ax.plot(ages,y_uncond[:Tmax,ik],marker=marker_list[0],color=color_list[0],label='unconditional (simulated)')
        ax.plot(ages,y_cond[:Tmax,ik],marker=marker_list[1],color=color_list[1],label='within group (weighted)')

        # title
        tit = ''
        if prefix != '': 
            tit = f'{prefix}_within_vs_between_k{k}_{mom_name}'
        ax.set_title(tit)

        # details
        ax.legend()
        ax.set_ylabel(mom_name + '$(y_{it+'+str(k)+'} - y_{it})$') 
        ax.set_xticks(np.arange(min(ages),max(ages),5))

        # save
        if prefix != '': 
            fig.tight_layout()
            save(fig,ax,tit)
            plt.close(fig)
        else: 
            plt.show()    