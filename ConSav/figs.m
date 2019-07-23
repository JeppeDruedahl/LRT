classdef figs
methods(Static)
    
%%%%%%%%%%%%%%%%%%%%%%%
% 1. policy functions %
%%%%%%%%%%%%%%%%%%%%%%%

function [] = polfuncs(vary,vary_latex,par,sol,ts,d,i_N,i_P)
        
    fig = figure('Name',sprintf('%s_t_%d_%d_d%d_N%d_P%d',...
        vary,ts(1),ts(end),d,i_N,i_P));
    hold('on');
    
    for t = ts
        x = sol.M_ast{d,t}(:,i_N,i_P); 
        y = sol.(vary){d,t}(:,i_N,i_P);
        h = plot(x,y,'o','MarkerSize',2,'DisplayName',sprintf('$t = %d$',t));
        set(h, 'MarkerFaceColor', get(h, 'Color'));       
    end

    % layout
    xlabel('$M_t$');    
    ylabel(vary_latex);
    legend('Location','best');
    box('on');
    grid on;
    
    funs.printfig(par,fig);
    
    close all;

end

function [] = all_polfuncs(par,sol)
    
    ts = [1 10 20 25 28 29];
    for d = unique(ceil(linspace(1,double(min(par.Nd(ts))),3)))
    for i_N = [1 par.NN]
    for i_P = unique([1 ceil(par.NP/2) par.NP])
        figs.polfuncs('C_ast','$C_t$',par,sol,ts,d,i_N,i_P);
        figs.polfuncs('V_ast','',par,sol,ts,d,i_N,i_P);        
    end
    end
    end
    
    close all;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%
% 2. simulation figure %
%%%%%%%%%%%%%%%%%%%%%%%%

function [] = d_hist(ts,par,sim) 
   
    for t = ts
        
        fig = figure('Name',sprintf('d_hist_t%d',t));

        histogram(sim.d(:,t),par.Nd(t))

        xlabel('group');
        ylabel('freq.');
        box('on');
        grid on;

        funs.printfig(par,fig);

    end
    
    close all;

end
    
function [] = lcp(vary,method,vary_latex,par,sim)
    
    if isnumeric(method)
        fig = figure('Name',sprintf('lcp_%s_p%d',vary,method));
    else
        fig = figure('Name',sprintf('lcp_%s_%s',vary,method));        
    end
    hold('on') ;

    x = (1:par.simT)+par.age_min-1;

    if strcmp(method,'mean') == 1
        y_sim  = mean(sim.(vary),1);  
    else
        y_sim  = prctile(sim.(vary),method,1);  
    end
    
    color = hex2rgb('1f77b4')/255;
    h = plot(x,y_sim,'-o','MarkerSize',5, 'Linewidth', 1.5, 'Color', color);
    set(h, 'MarkerFaceColor', get(h, 'Color'));

    if ~isnumeric(method)
        if contains('CY',vary)
            ylim([0 2])
        elseif strcmp(vary,'A')
            ylim([0 5])
        elseif strcmp(vary,'N')
            ylim([0 7])        
        end    
    end
    
    xlabel('age','FontSize',16);
    ylabel(vary_latex,'FontSize',16);
   
    ax = ancestor(h,'axes');
    Yaxis = ax.YAxis;
    Yaxis.FontSize = 16;
    Xaxis = ax.XAxis;
    Xaxis.FontSize = 16;
    
    box('on');
    grid on;
    
    funs.printfig(par,fig);
    
    close all;

end
function [] = lcp_compare(vary,method,vary_latex,par,sim,sim_alt,postfix,labels)
    
    if isnumeric(method)
        fig = figure('Name',sprintf('lcp_%s_%s_p%d',postfix,vary,method));
    else
        fig = figure('Name',sprintf('lcp_%s_%s_%s',postfix,vary,method));        
    end
    hold('on') ;

    x = (1:par.simT)+par.age_min-1;

    if strcmp(method,'mean') == 1
        y_sim  = mean(sim.(vary),1); 
        if ~strcmp(postfix,'data')
            y_sim_norisk  = mean(sim_alt.(vary),1);          
        end
    elseif strcmp(method,'var') == 1
        y_sim  = var(sim.(vary),0,1); 
        if ~strcmp(postfix,'data')
            y_sim_norisk  = var(sim_alt.(vary),1);          
        end        
    else
        y_sim  = prctile(sim.(vary),method,1);  
        if ~strcmp(postfix,'data')
            y_sim_norisk  = prctile(sim_alt.(vary),method,1);          
        end
    end
    if strcmp(postfix,'data')
        y_sim_norisk  = sim_alt;
    end
    
    if strcmp(postfix,'data')
        
        h = plot(x,y_sim_norisk,'-s','MarkerSize',5,'Linewidth',1.5,...
            'Color','black','DisplayName',labels{2});
        set(h, 'MarkerFaceColor', get(h, 'Color'));
        
        color = hex2rgb('1f77b4')/255;
        h = plot(x,y_sim,'-o','MarkerSize',5,'Linewidth',1.5,...
            'Color',color,'DisplayName',labels{1});
        set(h, 'MarkerFaceColor', get(h, 'Color'));
        
    else    
        
        color = hex2rgb('1f77b4')/255;
        h = plot(x,y_sim,'-o','MarkerSize',5,'Linewidth',1.5,...
            'Color',color,'DisplayName',labels{1});
        set(h, 'MarkerFaceColor', get(h, 'Color'));
        
        color = hex2rgb('ff7f0e')/255;
        h = plot(x,y_sim_norisk,'-o','MarkerSize',5,'Linewidth', 1.5,...
            'Color',color,'DisplayName',labels{2});
        set(h, 'MarkerFaceColor', get(h, 'Color'));
    
    end
    
    ax = ancestor(h,'axes');
    Yaxis = ax.YAxis;
    Yaxis.FontSize = 16;
    Xaxis = ax.XAxis;
    Xaxis.FontSize = 16;
    
    if ~isnumeric(method)
        if contains('CY',vary)
            ylim([0 2])
        elseif strcmp(vary,'A')
            ylim([0 3])
        elseif strcmp(vary,'N')
            ylim([0 7])        
        end    
    end
    xlabel('age','FontSize',20);
    ylabel(vary_latex,'FontSize',20);
    box('on');
    grid on;
    lgd = legend('Location','best');
    lgd.FontSize = 14;
    
    funs.printfig(par,fig);
    
    close all;

end

end
end