classdef funs
methods(Static)
    
function [] = layout()

    % set layout parameters
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter','latex');
    set(groot, 'defaultAxesFontSize', 12); 
    warning('off','all')

end  
function [x,fval,exitflag,output] = fmincon(goal,x0,lb,ub,options)
    
    [x,fval,exitflag,output] = fmincon(goal,x0,[],[],[],[],lb,ub,[],options);

end
function y = vec(x)
    y = x(:);
end
function [x, w] = GaussHermite(n)

    i   = 1:n-1;
    a   = sqrt(i/2);
    CM  = diag(a,1) + diag(a,-1);
    [V, L]   = eig(CM);
    [x, ind] = sort(diag(L));
    V       = V(:,ind)';
    w       = sqrt(pi) * V(:,1).^2;

end
function [x, w] = GaussHermite_lognorm(sigma,n,correct_mean)

    [x,w] = funs.GaussHermite(n);

    if correct_mean
        x = exp(x*sqrt(2)*sigma-0.5*sigma^2);
    else
        x = exp(x*sqrt(2)*sigma);
    end
    w = w./sqrt(pi);

    % assert a mean of one
    if correct_mean
        assert(1-sum(w.*x) < 1e-8)
    end
    
end    
function x = nonlinspace(lo,hi,n,phi)
    % recursively constructs an unequally spaced grid.
    % phi > 1 -> more mass at the lower end of the grid.
    % lo can be a vector (x then becomes a matrix).

    x      = NaN(n,length(lo));
    x(1,:) = lo;
    for i = 2:n
        x(i,:) = x(i-1,:) + (hi-x(i-1,:))./((n-i+1)^phi);
    end

end
function [] = printfig(par,figin)

    fig = figure(figin);
    fig.PaperUnits = 'centimeters';   
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0 0 16 12];
    fig.PaperSize = [16 12];

    folder = ['figs_tabs\' par.prefix];
    if exist(folder,'dir') == 0
        mkdir(folder);
    end
    filename = [folder '\' get(fig,'name') ''];
    print('-dpdf',['' filename '.pdf']);

end
function [] = printtable(rownames,table,filename)
    
    file = fopen(sprintf('%s/%s',pwd,filename),'w');
    rows = size(table,1);
    cols = size(table,2);
    for row = 1:rows     
        fprintf(file,'%s &',rownames{row});
    for col = 1:cols
        if numel(table{row,col}) > 0
            fprintf(file,'%s',table{row,col});
        end
        if col < cols
            fprintf(file,' & ');
        else
            fprintf(file,' \\\\ \n');
        end
    end
    end
    fclose(file);
    
end


end
end